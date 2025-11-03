import yaml
import logging
import os
import datetime
import subprocess
import shutil
import sys
import itertools
import copy

from define import Ingredient, Constraint, Role, Result, get_constraints_idx, reduce_guests, print_reduced
from dock import dock
from utils import read_xyz, merge_xyz, split_docker_results, split_multi_pdb, split_multi_xyz, write_multi_pdb, append_scores, get_atom_count
from evaluate import filter_conformations
from arrange import arrange_guests
from optimise import optimise
from drOrca import make_orca_input

########################
## LOGGING
########################

def setup_logging(workdir):
    # Remove existing handlers
    root_logger = logging.getLogger()
    for handler in list(root_logger.handlers):
        root_logger.removeHandler(handler)

    # Log to both console and file
    handlers = [logging.StreamHandler()]  # console
    handlers.append(logging.FileHandler(os.path.join(workdir, "log.txt")))  # file

    logging.basicConfig(
        format="%(asctime)s - %(levelname)s - %(message)s",
        level=logging.DEBUG,
        handlers=handlers,
    )

########################
## SETUP
########################

def setup_config(configPath):
    if not os.path.isfile(configPath):
        sys.exit(f"ERROR: File does not exist at {configPath}")
    # Load the YAML file
    with open(configPath, "r") as file:
        config = yaml.safe_load(file)
    logging.info(f"Loaded configuration from {configPath}\n")
    misc = config.get("misc", {})
    workdir = misc.get("workdir", ".")
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    outdir = os.path.join(workdir, f"output_{timestamp}")
    os.makedirs(outdir, exist_ok=True)
    qmMethod = misc.get("qmMethod", "XTB2")
    quick = misc.get("quick", False)
    nprocs = misc.get("nprocs", 8)
    orcapth = misc.get("orcapth", "orca")
    return config, outdir, qmMethod, quick, nprocs, orcapth

def setup_ingredients(config):
    ingredients_cfg = config.get("ingredients", [])
    roles_cfg = config.get("roles", [])

    logging.debug(f"""Ingredients are:
                {ingredients_cfg}\n""")
    logging.debug(f"""Roles are:
                {roles_cfg}\n""")

    ingredient_map = {}  # Map ingredient names to objects for role processing
    for ing in ingredients_cfg:
        roles_dict = ing.get('roles', {})
        ingredient_obj = Ingredient(
            path=ing['path'],
            eopt=None,
            einter=None,
            charge=ing['charge'],
            multiplicity=ing['multiplicity'],
            roles=roles_dict,
            name=ing['name']
        )
        ingredient_map[ing['name']] = ingredient_obj
    logging.debug(f"""Ingredient map is:
                {ingredient_map}\n""")

    # Create Role objects from YAML
    role_objects = []
    for role in roles_cfg:
        # Map candidates to Ingredient objects
        host_candidates = []
        for host_cand in role["host_candidates"]:
            if "guests_of" in host_cand:
                host_candidates = [ingredient_map[cand['name']] for cand in host_cand["guests_of"]["guest_candidates"]]
            else:
                host_candidates = [ingredient_map[cand['name']] for cand in role['host_candidates']]
        guest_candidates = [ingredient_map[cand['name']] for cand in role['guest_candidates']]
        
        # Handle constraints if present
        constraints_data = role.get('constraints', [])
        constraints = []
        for cons in constraints_data:
            constraint = Constraint(
                guestIdx=cons['guestIdx'],
                guestType=cons.get('guestType'),
                hostIdx=cons['hostIdx'],
                hostType=cons.get('hostType'),
                val=cons['val'],
                force=cons['force']
            )
            constraints.append(constraint)
        
        # Create Role object
        role_obj = Role(
            title=role['name'],
            priority=role['priority'],
            guests=guest_candidates,
            host=host_candidates,
            constraints=constraints
        )
        role_objects.append(role_obj) 
    logging.debug(f"""Role objects are:
                {role_objects}""")

    return role_objects, ingredient_map

def process_reindexing(protonated_structures, reindex_reference_path, reindex_output_dir, ing_name, ing_id, logger):
    successful = []
    failed = []
    for idx, path in enumerate(protonated_structures):
        reindexed_pdb_path = reindex(reindex_reference_path, path, reindex_output_dir, logger=logger)
        if reindexed_pdb_path is None:
            if logger:
                logging.info(f"Reindexing failed for pose {idx} of ingredient {ing_name}")
            continue
        successful.append(reindexed_pdb_path)

    # Check if all poses failed
    if not successful:
        failed.append(ing_id)
        if logger:
            logging.warning(f"No poses were successfully reindexed for ingredient {ing_name}")
    else:
        if logger:
            logging.info(f"For ingredient {ing_name}, {len(successful)} out of {len(protonated_structures)} poses were successfully reindexed")

    return successful, failed

########################
## HELPER FUNCTIONS
########################

def calculate_charge(charges):
    return int(sum(charges))
    
def calculate_multiplicity(multiplicities):
    # Assume molecules are weakly or non-interacting
    spin = int(sum(m - 1 for m in multiplicities) / 2)
    multiplicity = 2*spin + 1
    return multiplicity

def write_docking_input(role_name, guest, host, constraints, workdir, qmMethod, quick, nprocs):
    inp_file_path = os.path.join(workdir, f"dock")
    title = f"Docking {'(quick)' if quick else ''}\n # role: {role_name}\n # host: {host.name}\n # guest: {guest.name}"

    hostPathPDB = host.path
    hostPathXYZ = hostPathPDB.replace(".pdb", ".xyz")
    # TODO use openbabel to convert PDB to XYZ whether XYZ exists or not - to ensure atom index consistency
    # TODO need to install obabel first
    if not os.path.isfile(hostPathXYZ):
        logging.error(f"""Host XYZ file not found at {hostPathXYZ}.
                      In the future we will implement autmoatic obabel conversion.
                      For now, please convert the PDB to XYZ to ensure atom index consistency.""")
        
    guestPathPDB = guest.path # TODO this needs to work with host being product of previous docking
    guestPathXYZ = guestPathPDB.replace(".pdb", ".xyz")
    if not os.path.isfile(guestPathXYZ): # TODO
        logging.error(f"""Guest XYZ file not found at {guestPathXYZ}.
                      In the future we will implement autmoatic obabel conversion.
                      For now, please convert the PDB to XYZ to ensure atom index consistency.""")
    
    simpleInputLine = [qmMethod]
    if quick:
        simpleInputLine.append("QUICKDOCK")
    inputFormat = "xyzfile"

    # TODO this needs to work with host being product of previous docking
    moleculeInfo = {"charge": calculate_charge([guest.charge, host.charge]),
                    "multiplicity": calculate_multiplicity([guest.multiplicity, host.multiplicity])}

    # "Bond bias potential": 
    # Atom number for GUEST FIRST, then HOST SECOND,
    # absolute indexing for each system independently
    biases = []
    for constraint in constraints:
        bias = {"atoms": [constraint.guestIdx, constraint.hostIdx], "val": constraint.val, "force": constraint.force}
        biases.append(bias)
    docker = {"guestPath": guestPathXYZ, "fixHost": True, "bias": biases}

    # Use the updated XYZ file for optimization
    make_orca_input(orcaInput=inp_file_path,
                    title=title,
                    simpleInputLine=simpleInputLine,
                    inputFormat=inputFormat,
                    inputFile=hostPathXYZ,
                    moleculeInfo=moleculeInfo,
                    parallelize=nprocs,
                    qmmm=None,
                    geom=None,
                    neb=None,
                    scf=None,
                    docker=docker)
    
    # Metadata returned to facilitate further docking, reusing products of earlier docking as hosts
    return inp_file_path, moleculeInfo['charge'], moleculeInfo['multiplicity'], [host, guest]

def run_docking(input, orcapath):
    output_file = f"{input}.out"
    with open(output_file, "w") as f:
        subprocess.run([orcapath, f"{input}.inp"], check=True, stdout=f, stderr=subprocess.STDOUT)

def process_docking_output(inp_file_path, curr_charge, curr_multiplicity, curr_ingredients, role_name, ingredient_map, logger):
    all_results = split_docker_results(f"{inp_file_path}.docker.struc1.all.optimized.xyz", logger)
    print(curr_ingredients[0].name, curr_ingredients[0].n_atoms)

    # TODO reuse previous Result objects... eopt = new_eopt, einter = einter + new_einter
    # TODO next docking should use all poses that satisfy constraints
    for result in all_results:
        (path, eopt, einter, df) = result
        # TODO evaluate output for satisfying constraints
        new_path = os.path.join(os.path.dirname(path), "passed", os.path.basename(path))
        os.makedirs(os.path.dirname(new_path), exist_ok=True)
        shutil.move(path, new_path)
        # Edit df to assign ingredient names to atoms
        df[["ING", "DISH"]] = None, None
        first_atom = 0
        # Don't want to be rewriting previous ingredients when reusing result of previous docking - this is just initial setup?
        for ing in curr_ingredients:
            last_atom = first_atom + ing.n_atoms
            df.loc[first_atom:last_atom-1, "ING"] = ing.name
            df.loc[first_atom:last_atom-1, "DISH"] = role_name
            first_atom = last_atom
        product_obj = Result(
            path=new_path,
            eopt=eopt,
            einter=einter,
            charge=curr_charge,
            multiplicity=curr_multiplicity,
            df=df
        )
        print(product_obj.path)
        print(product_obj.eopt)
        print(product_obj.einter)
        print(product_obj.charge)
        print(product_obj.multiplicity)
        print(product_obj.df)
        ingredient_map[role_name] = product_obj
        exit()

def process_docking(ing, host, outdir, dock_params, redocking=False):
    # Handles the docking, per-conformation protonation, reindexing, and deprotonation for a single ingredient.
    logging.info(f"Processing docking for {ing.name}{' (redocking)' if redocking else ''}...")
    # Create a temporary directory for per-conformation processing
    temp_processing_dir = os.path.join(outdir, f"temp_processing_{ing.name}")
    prot_output_dir = os.path.join(temp_processing_dir, f"protonation")
    reindex_output_dir = os.path.join(temp_processing_dir, f"reindexing")
    deprot_output_dir = os.path.join(temp_processing_dir, f"deprotonation")
    os.makedirs(prot_output_dir, exist_ok=True)
    os.makedirs(reindex_output_dir, exist_ok=True)
    os.makedirs(deprot_output_dir, exist_ok=True)

    reference_atom_count = get_atom_count(ing.path, logger) # Atom count of original guest molecule

    # Perform initial docking
    docked_path = dock(host, ing, outdir, dock_params, redocking, logger)
    if not os.path.exists(docked_path) or os.path.getsize(docked_path) == 0:
        logging.warning(f"Docked output for {ing.name} not found or is empty at {docked_path}. Skipping further processing.")
        return None
    
    # get rid of protonateion/deprotonation
    prot_path = docked_path

    # Prepare reference for reindexing (PDB version of original ingredient)
    reindex_reference_path = f"{os.path.splitext(ing.path)[0]}_ref.pdb"
    if not os.path.exists(reindex_reference_path):
        logging.warning("Reference ingredient in PDB format not found - using XYZ instead; consider adding a PDB reference with same atom indices and path as XYZ.")
        reindex_reference_path = ing.path

    # Split the protonated multi-PDB into individual PDB files
    protonated_structures = split_multi_pdb(prot_path, prot_output_dir, logger)
    # Reindex protonated PDB files individually
    successful, failed = process_reindexing(protonated_structures, reindex_reference_path, reindex_output_dir, ing.name, ing.id, logger)
    if failed:
        logging.error(f"Reindexing failed for {ing.name}. Investigate output manually.")
        logging.warning(f"Failed ingredients: {', '.join(failed)}")
        return None, failed
    # Merge protonated and reindexed PDB files into one multi-PDB file
    multi_prot_path = os.path.join(outdir, ing.name, f"{ing.name}_docked_prot.pdb")
    write_multi_pdb(successful, multi_prot_path)

    # Deprotonate specific atoms
    deprot_path, deprot_atom_count = deprotonate_selected(multi_prot_path, getattr(ing.indices, "deprotonate", []))
    if not os.path.exists(deprot_path):
        logging.error(f"Deprotonated file not found at expected path: {deprot_path}.")
    if deprot_atom_count != reference_atom_count:
        failed.append(ing.name)
        logging.error(f"Deprotonation failed for {ing.name}. Investigate output manually.")
        logging.warning(f"Failed ingredients: {', '.join(failed)}")
        return None, failed
    
    # Convert protonated, reindexed, deprotonated multi-PDB file to multi-XYZ file
    prot_reidx_deprot_pdb_path = deprot_path.replace(".pdb", ".xyz")
    cmd = [
        "obabel",
        "-ipdb", deprot_path,
        "-O", prot_reidx_deprot_pdb_path
        ]
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    # Merge this multi-XYZ guest file with the host structure
    final_merged_host_guest_path = os.path.join(outdir, f"{ing.name}.xyz")
    merge_xyz(host.path, prot_reidx_deprot_pdb_path, final_merged_host_guest_path)
    logging.info(f"Final merged host-guest file saved to: {final_merged_host_guest_path}")

    # Append scores from gnina docking
    scores_file_from_dock = os.path.join(os.path.dirname(docked_path), "scores.txt")
    if os.path.exists(scores_file_from_dock):
        append_scores(final_merged_host_guest_path, scores_file_from_dock, logger)
        logging.info(f"Appended scores to {final_merged_host_guest_path}")
    else:
        logging.warning(f"Scores file not found at {scores_file_from_dock}. Skipping score appending.")

    # Clean up temporary processing directory
    shutil.rmtree(temp_processing_dir)
    logging.debug(f"Cleaned up temporary directory: {temp_processing_dir}")

    return final_merged_host_guest_path, failed

def process_constraints(docked_path, ing, host, outdir, evaluate_backbone, dock_params, redock_params, logger):
    # Filter conformations that don't satisfy constraints - only valid written to filtered_path
    valid_structures, filtered_path = filter_conformations(docked_path, host.path, ing.name, ing.role_title, ing.indices, ing.constraints, evaluate_backbone, logger)

    # Repeat docking if no valid conformations found
    if len(valid_structures) == 0:
        logging.warning(f"No poses that satisfy constraints found on first docking attempt. Redocking with {redock_params}...")
        # Update dock_params with redocking parameters
        redock_dock_params = dock_params.copy()
        redock_dock_params.update(redock_params)
        # Repeat docking with redocking parameters
        redocked_path, failed = process_docking(ing, host, outdir, redock_dock_params, redocking=True)

        # Repeat filtering conformations that don't satisfy constraints
        valid_structures, filtered_path = filter_conformations(redocked_path, host.path, ing.name, ing.role_title, ing.indices, ing.constraints, evaluate_backbone, logger)
        if len(valid_structures) == 0:
            logging.error(f"No poses that satisfy constraints {ing.constraints} for {ing.name}, role: {ing.role_title} found on repeat docking.")
            return
    ing.update_conformations(filtered_path)
    logging.info(f"Filtered conformations are saved in {filtered_path}")
    return filtered_path

def expand_constraints(guest, host, constraint):
    # Filter tuples whose activity matches the requested activity name
    guest_matches = [(i, row["ATOM_NAME"], row["ROLE"]) for i, row in guest.df.iterrows() if constraint.guestIdx in row["ROLE"]]
    host_matches = [(i, row["ATOM_NAME"], row["ROLE"]) for i, row in host.df.iterrows() if constraint.hostIdx in row["ROLE"]]

    keep = []
    for g_item, h_item in itertools.product(guest_matches, host_matches):
        # g_item and h_item are tuples like (index, name, activity)
        g_idx = g_item[0]  # numeric index in the PDB/molecule
        h_idx = h_item[0]
        keep.append((g_idx, h_idx, constraint.val, constraint.force))
    return keep

def expand_role_combinations(role):
    """
    Expand roles with multiple hosts/guests into all (host, guest) pairs.
    For each host/guest pair, expand constraints so that for every possible
    combination of concrete atom-pairs (one per constraint) we produce a unique Role.
    """
    expanded_roles = {}

    # Loop all candidate host/guest pairs
    for host, guest in itertools.product(role.host, role.guests):
        # For each constraint compute the list of concrete options (g_idx, h_idx, val, force)
        options_per_constraint = []
        for constraint in role.constraints or []:
            opts = expand_constraints(guest, host, constraint)
            # If a constraint has no matching atoms, this role/host/guest pair can't satisfy it:
            if not opts:
                # skip this host/guest pair entirely
                options_per_constraint = []
                break
            options_per_constraint.append(opts)

        # If there were no constraints, create a single entry (no-constraint case)
        if not options_per_constraint and (not role.constraints):
            # Create single role copy with no concrete constraints
            new_role = copy.deepcopy(role)
            new_role.host = host
            new_role.guests = guest
            new_role.constraints = []
            role_key = f"{role.title}_{host.name}_{guest.name}_no_constraints"
            expanded_roles[role_key] = new_role
            continue

        # If any constraint had zero options, skip (cannot satisfy constraints)
        if not options_per_constraint:
            continue

        # Cartesian product across lists of options (one option per constraint)
        # If one ingredient has multiple atoms of same flavour (e.g. multiple h_donor),
        # this creates multiple combinations, one for each of the atoms
        for combo_idx, combo in enumerate(itertools.product(*options_per_constraint)):
            # combo is a tuple of option tuples, one per constraint
            # Build new constraint objects (or simple tuples) representing these concrete constraints
            concrete_constraints = []
            for orig_constraint, option in zip(role.constraints, combo):
                g_idx, h_idx, val, force = option
                # Create a copy of the original constraint but with resolved atom indices
                new_cons = copy.deepcopy(orig_constraint)
                # Overwrite guestIdx/hostIdx with the concrete atom indices
                new_cons.guestIdx = g_idx
                new_cons.hostIdx = h_idx
                # Ensure val/force are set to the chosen option (val likely same as orig)
                new_cons.val = val
                new_cons.force = force
                concrete_constraints.append(new_cons)

            # Make new role copy and set host/guest to the concrete candidates and constraints
            new_role = copy.deepcopy(role)
            new_role.host = host
            new_role.guests = guest  # TODO rename to "guest" considering the 1-to-1 expansion
            new_role.constraints = concrete_constraints

            # Unique key: include the combination index so duplicates don't clobber
            role_key = (str(role.title), str(host.name), str(guest.name), str(combo_idx)) # is a tuple
            expanded_roles[role_key] = new_role

    return expanded_roles

########################
## MAIN LOOP
########################

def main(args):
    # Read config file
    if not args.config:
        args.config = "/home/mchrnwsk/prometheozyme/config.yaml"
    config, outdir, qmMethod, quick, nprocs, orcapath = setup_config(args.config)
    setup_logging(outdir)

    # Prepare ingredients and roles from the configuration    
    roles, ingredient_map = setup_ingredients(config)

    # Run ORCA DOCKER for each role, exhaustive for host/guest mix
    failed = []
    for role in roles:
        expanded_roles = expand_role_combinations(role)
        logging.debug(f"""Expanded ingredient and constraint combinations for
                      role title: {role.title}
                      number of combinations: {len(expanded_roles)}
                      combination keys: {list(expanded_roles.keys())}""")
        for role_key, role_desc in expanded_roles.items():
            logging.info(f"Processing role: {role_key}")

            # Output is organised as role/host/guest/constraints_combination
            workdir = os.path.join(outdir, *role_key)
            os.makedirs(workdir, exist_ok=True)
            inp_file_path, curr_charge, curr_multiplicity, curr_ingredients = write_docking_input(role_key[0], role_desc.guests, role_desc.host, role_desc.constraints, workdir, qmMethod, quick, nprocs)

            # run_docking(inp_file_path, orcapath)
            inp_file_path = "runs/output_2025-11-02_21-36-01/base1/sub/his/0/dock"
            process_docking_output(inp_file_path, curr_charge, curr_multiplicity, curr_ingredients, role_key[0], ingredient_map, logging)



    sys.exit(1)
        #if ing.name == 'host':
        #    continue
        #logging.info(f"Running gnina for ingredient: {ing.name}")
        #outpath, failed_i = process_docking(ing, host, outdir, dock_params)
        #if failed_i:
        #    failed.append(failed_i)
    logging.info(f"Docking completed. Results saved in {outdir}\n")
    
    # Evaluate constraints for each ingredient (list of unique guests)
    # for ingredient in unique_guests_constraints, delete in copy conformations that don't satisfy constraints
    for ing in unique_guests_constraints:
        logging.info(f"\nEvaluating constraints for ingredient: {ing.name}, role: {ing.role_title}")

        docked_path = os.path.join(outdir, f"{ing.name}.xyz")
        if not os.path.exists(docked_path):
            logging.error(f"Docked output for {ing.name} not found at {docked_path}. Skipping constraint evaluation.")
            continue

        # Filter conformations that don't satisfy constraints
        process_constraints(docked_path, ing, host, outdir, evaluate_backbone, dock_params, redock_params, logger)
        
    # Satisfy roles according to their priorities
    sorted_final_arrangements = arrange_guests(roles, unique_guests_constraints, host_atom_count, host_coords, host_atom_types, outdir, logger)

    for arr in sorted_final_arrangements:
        failed_ids = [guest["obj"] for guest in arr["guests_info"] if guest["obj"] in failed]

        if failed_ids:
            logging.warning(f"Arrangement at path {arr['path']} includes failed guests with IDs: {failed_ids}. Protonate/reindex/deprotonate manually before proceeding.")
        else:
            logging.debug(f"Arrangement at path {arr['path']} has no failed guests.")
        # simple constrained optimisation
        arr_optimised, arr_reactant = optimise(arr, host_atom_count, ingredient_map, backbone, logger)
    #    # prepare product
    #    arr_product = prepare_product(arr, arr_reactant)

    #for arrangement in final_arrangements_list:
        # run NEB

    logging.info("""
          Main script complete!
          ⚝⭒٭⋆⚝⭒٭⋆⚝⭒٭⋆⚝⭒٭⋆⚝⭒٭⋆⚝⭒٭⋆
          """)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Process some ingredients and roles.")
    parser.add_argument('--config', type=str, default=None, help="Path to the configuration YAML file")
    args = parser.parse_args()
    main(args)