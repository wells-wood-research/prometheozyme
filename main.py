import yaml
import logging
import os
import datetime
import subprocess
import shutil
import sys
import itertools
import copy
import string

from define import Ingredient, Constraint, Role, col_order, col_types
from utils import split_docker_results
from evaluate import filter_conformations
from drOrca import make_orca_input
import pdbUtils

########################
## SETUP
########################

def setup(configPath):
    config, outdir, orca = setup_config(configPath)
    setup_logging(outdir)
    # Prepare ingredients and roles from the configuration    
    courses, ingredients = setup_ingredients(config) # TODO renamed roles --> courses, ingredient_map --> ingredients here
    return config, outdir, orca, courses, ingredients

def setup_config(configPath):
    if not os.path.isfile(configPath):
        sys.exit(f"ERROR: File does not exist at {configPath}")
    # Load the YAML file
    with open(configPath, "r") as file:
        config = yaml.safe_load(file)
    logging.info(f"Loaded configuration from {configPath}\n")

    # Get miscellaneous parameters
    misc = config.get("misc", {})
    workdir = misc.get("workdir", ".")

    # Get orca docker parameters
    orca = config.get("orca", {})

    # Setup output dir
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    outdir = os.path.join(workdir, f"output_{timestamp}")
    os.makedirs(outdir, exist_ok=True)

    return config, outdir, orca

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

def setup_ingredients(config):
    ingredients_cfg = config.get("ingredients", [])
    courses_cfg = config.get("courses", [])

    logging.debug(f"""Ingredients are:
                {ingredients_cfg}\n""")
    logging.debug(f"""Courses are:
                {courses_cfg}\n""")

    # Create 'Ingredient' objects
    # TODO is ingredients map really useful?
    # TODO are all fields useful?
    ingredients = {}
    for ing in ingredients_cfg:
        ingredient_obj = Ingredient(
            name=ing['name'],
            path=ing['path'],
            eopt=0,
            einter=0,
            charge=ing['charge'],
            multiplicity=ing['multiplicity'],
            flavours=ing.get('flavours', {}),
            df=None
        )
        ingredients[ing['name']] = ingredient_obj
    logging.debug(f"""Ingredient map is:
                {ingredients}\n""")

    # Create 'Course' objects
    course_objects = {}
    for course in courses_cfg:
        if course['name'] == 'init':
            host_candidates = None
            guest_candidates = [ingredients["sub"]] # list for consistent typing
            constraints = None
        else:
            # Map candidates to Ingredient objects
            host_candidates = course['host_candidates']['name']
            guest_candidates = [ingredients[cand['name']] for cand in course['guest_candidates']]
            
            # Handle constraints if present
            constraints_data = course.get('constraints', [])
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
        
        # Create Course object
        course_obj = Course(
            name=course['name'],
            host=host_candidates,
            guests=guest_candidates,
            constraints=constraints
        )
        course_objects[course_obj.title] = course_obj 
    logging.debug(f"""course objects are:
                {course_objects}""")

    return course_objects, ingredients

########################
## HELPER FUNCTIONS
########################

def assign_chain_ids(df):
    """
    Ensures that every unique (ING, DISH) pair has a unique CHAIN_ID.
    Keeps existing CHAIN_IDs if consistent, assigns new ones in A, B, C... order otherwise.
    """
    # Get all unique (ING, DISH) pairs in the order they appear
    unique_pairs = []
    for _, row in df.iterrows():
        pair = (row["ING"], row["DISH"])
        if pair not in unique_pairs:
            unique_pairs.append(pair)
    
    # Generate chain labels A, B, C, ..., AA, AB, ...
    def get_chain_label(i):
        letters = string.ascii_uppercase
        label = ""
        while True:
            i, r = divmod(i, 26)
            label = letters[r] + label
            if i == 0:
                break
            i -= 1
        return label

    # Assign chain labels to each unique (ING, DISH)
    chain_map = {pair: get_chain_label(i) for i, pair in enumerate(unique_pairs)}
    
    # Apply mapping to the DataFrame
    df["CHAIN_ID"] = [chain_map[(ing, dish)] for ing, dish in zip(df["ING"], df["DISH"])]
    
    return df

def calculate_charge(charges):
    return int(sum(charges))
    
def calculate_multiplicity(multiplicities):
    # Assume molecules are weakly or non-interacting
    spin = int(sum(m - 1 for m in multiplicities) / 2)
    multiplicity = 2*spin + 1
    return multiplicity

def dock(outdir, role_key, role_desc, orca):
    # Read orca parameters
    orcapath = orca.get("orcapth", "./orca")
    qmMethod = orca.get("qmMethod", "XTB2")
    strategy = orca.get("strategy", "normal")
    optLevel = orca.get("optLevel", "sloppyopt")
    nOpt = orca.get("nOpt", 5)
    gridExtent = orca.get("gridExtent", 15)
    nprocs = orca.get("nprocs", 8) 

    # Output of each docking step is organised as role/guest/constraints_combination
    workdir = os.path.join(outdir, *role_key)
    os.makedirs(workdir, exist_ok=True)

    inp_file_path, curr_charge, curr_multiplicity = write_docking_input(role_key[0], role_desc.guests, role_desc.host, role_desc.constraints, workdir, qmMethod, strategy, optLevel, nOpt, gridExtent, nprocs)

    run_docking(inp_file_path, orcapath)
    results_map = process_docking_output(inp_file_path, curr_charge, curr_multiplicity, role_desc.guests, role_desc.host, role_key, logging)
    return results_map

def write_docking_input(role_name, guest, host, constraints, workdir, qmMethod, strategy, optLevel, nOpt, gridExtent, nprocs):
    inp_file_path = os.path.join(workdir, f"dock")
    title = f"ORCA DOCKER: Automated Docking Algorithm for\n # role: {role_name}\n # host: {host.name}\n # guest: {guest.name}"

    # Prepare host and guest molecular structure files as XYZ
    hostPathPDB = host.path
    hostPathXYZ = hostPathPDB.replace(".pdb", ".xyz")
    # TODO convert PDB to XYZ whether XYZ exists or not - to ensure atom index consistency
    if not os.path.isfile(hostPathXYZ):
        logging.error(f"""Host XYZ file not found at {hostPathXYZ}.
                      In the future we will implement autmoatic obabel conversion.
                      For now, please convert the PDB to XYZ to ensure atom index consistency.""")
    guestPathPDB = guest.path
    guestPathXYZ = guestPathPDB.replace(".pdb", ".xyz")
    if not os.path.isfile(guestPathXYZ):
        logging.error(f"""Guest XYZ file not found at {guestPathXYZ}.
                      In the future we will implement autmoatic obabel conversion.
                      For now, please convert the PDB to XYZ to ensure atom index consistency.""")
    
    simpleInputLine = [qmMethod]
    inputFormat = "xyzfile"

    # charge and multiplicity of host only - guest is defined under %DOCKER GUESTCHARGE and GUESTMULT
    moleculeInfo = {"charge": host.charge, "multiplicity": host.multiplicity}

    # Define bond bias potential to impose restraints of distance between guest and host atoms
    # Must be defined by atom numbers for GUEST FIRST, then HOST SECOND, absolute indexing for each molecule independently
    # For example, to add a bond bias between atom 2 from the GUEST and atom 19 from the HOST, write: BIAS { B 2 19 } END 
    biases = []
    for constraint in constraints:
        bias = {"atoms": [constraint.guestIdx, constraint.hostIdx], "val": constraint.val, "force": constraint.force}
        biases.append(bias)
    
    docker = {"guestPath": guestPathXYZ, "guestCharge": guest.charge, "guestMultiplicity": guest.multiplicity, "fixHost": True, "bias": biases, "strategy": strategy, "optLevel": optLevel, "nOpt": nOpt, "gridExtent": gridExtent}

    # Use the updated XYZ file for optimization
    make_orca_input(orcaInput=inp_file_path,
                    title=title,
                    simpleInputLine=simpleInputLine,
                    inputFormat=inputFormat,
                    inputFile=hostPathXYZ,
                    moleculeInfo=moleculeInfo,
                    parallelize=nprocs,
                    docker=docker)
    
    # Metadata returned to facilitate further docking, reusing products of earlier docking as hosts
    return inp_file_path, calculate_charge([guest.charge, host.charge]), calculate_multiplicity([guest.multiplicity, host.multiplicity])

def run_docking(input, orcapath):
    output_file = f"{input}.out"
    with open(output_file, "w") as f:
        subprocess.run([orcapath, f"{input}.inp"], check=True, stdout=f, stderr=subprocess.STDOUT)

def process_docking_output(inp_file_path, curr_charge, curr_multiplicity, guest, host, role_key, logger):
    role_name = role_key[0]
    all_results = split_docker_results(f"{inp_file_path}.docker.struc1.all.optimized.xyz", logger)
    results_map = {}

    # TODO reuse previous Result objects... eopt = new_eopt, einter = einter + new_einter
    # TODO next docking should use all poses that satisfy constraints
    c = 0
    for result in all_results:
        (path, eopt, einter, df) = result
        # TODO evaluate output for satisfying constraints
        # TODO after optmisation there might be duplicated results - need to remove
        new_path = os.path.join(os.path.dirname(path), f"struct{c}", "host.xyz")
        os.makedirs(os.path.dirname(new_path), exist_ok=True)
        shutil.move(path, new_path)
        # Edit df to assign ingredient names to atoms
        df[["ROLE", "ING", "DISH"]] = None, None, None
        n_host, n_guest = host.n_atoms, guest.n_atoms
        # transfer original labels

        # TODO probably host df here is just the constrained element, not all
        
        for col in ["ATOM", "ATOM_ID", "ATOM_NAME", "RES_NAME", "CHAIN_ID", "RES_ID", "OCCUPANCY", "BETAFACTOR", "ELEMENT", "ROLE", "ING"]:
            df.loc[:n_host - 1, col] = host.df[col].values
            df.loc[n_host:n_host + n_guest - 1, col] = guest.df[col].values
        # transfer dish labels
        df.loc[:n_host - 1, "DISH"] = host.df["DISH"].values
        df.loc[n_host:n_host + n_guest - 1, "DISH"] = role_name

        # Make sure each molecule gets a new CHAIN_ID so that ATOM_IDs can be individually 1-based per each molecule
        print(df)
        df = assign_chain_ids(df)

        # TODO ATOM_ID should probably continue... of chain id should be changed...
        df = df[col_order].astype(col_types)

        # Save PDB
        df_to_save = df.loc[:, ~df.columns.isin(["ROLE", "ING", "DISH"])]
        pdbUtils.df2pdb(df_to_save, new_path.replace(".xyz", ".pdb"))
        print(f"Saved to {new_path.replace('.xyz', '.pdb')}")
        
        product_obj = Ingredient(
            path=new_path,
            eopt=eopt,
            einter=einter,
            charge=curr_charge,
            multiplicity=curr_multiplicity,
            df=df
        )
        new_role_key = role_key + (str(c), )
        results_map[new_role_key] = product_obj
        c+=1
    return results_map

def expand_constraints(guest, host, constraint):
    # Filter tuples whose activity matches the requested activity name
    print("For guest:")
    print(guest.df)
    print(constraint.guestIdx)
    print()
    guest_matches = [(i, row["ATOM_NAME"], row["ROLE"]) for i, row in guest.df.iterrows() if constraint.guestIdx in row["ROLE"]]
    print("For host:")
    print(host.df)
    print(constraint.hostIdx)
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
            role_key = (str(role.title), str(guest.name), str(combo_idx)) # is a tuple
            expanded_roles[role_key] = new_role

    return expanded_roles

########################
## MAIN LOOP
########################

def main(args):
    # Read config file
    if not args.config:
        # TODO del when finished testing
        args.config = "/home/mchrnwsk/prometheozyme/config.yaml"
    outdir, orca, roles, ingredients = setup(args.config)

    # Prepare the substrate as a fake result of first docking
    prev_results = []
    init_host = roles["init"].guests[0]
    init_host_dir = os.path.join(outdir, "init")
    init_host_path = os.path.join(init_host_dir, "host.pdb") # TODO do you want "prev_results" to be PDB or XYZ based? DO you want to always (!) write both XYZ and PDB?
    os.makedirs(init_host_dir)
    shutil.copy(init_host.path, init_host_path)
    # TODO assumes there's xyz present - need to expect user only to prepare PDBs
    shutil.copy(init_host.path.replace(".pdb", ".xyz"), init_host_path.replace(".pdb", ".xyz"))
    init_host.path = init_host_path
    prev_results.append(init_host)

    # Run ORCA DOCKER for each role, exhaustive for host/guest mix
    for i, (role_name, role) in enumerate(roles.items()):
        if role_name == "init":
            continue
        for result in prev_results:
            new_host_ing_for_docking = result
            new_host_df = new_host_ing_for_docking.df
            new_host_df_for_constraints = new_host_df[new_host_df["DISH"] == role.host]
            new_host_name = new_host_df_for_constraints["ING"].unique()[0]
            new_host_ing_for_constraints = ingredients[new_host_name]
            new_host_ing_for_constraints.df = new_host_df_for_constraints
            new_host_ing_for_constraints.path = result.path

            new_role = copy.deepcopy(role)
            new_role.host = [new_host_ing_for_constraints]

            expanded_roles = expand_role_combinations(new_role)
            logging.debug(f"""Expanded ingredient and constraint combinations for
                        role title: {role.title}
                        number of combinations: {len(expanded_roles)}
                        combination keys: {list(expanded_roles.keys())}""")

            new_prev_results = []
            for role_key, role_desc in expanded_roles.items():
                logging.info(f"Processing role: {role_key}")
                role_desc.host = new_host_ing_for_docking
                outdir = os.path.dirname(result.path)
                results_map = dock(outdir, role_key, role_desc, orca)
                [new_prev_results.append(results_map[key]) for key in results_map.keys() if key[:3] == role_key]
        prev_results = new_prev_results

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