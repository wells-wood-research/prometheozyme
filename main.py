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
            name=ing['name'],
            df=None
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
    return inp_file_path, moleculeInfo['charge'], moleculeInfo['multiplicity']

def run_docking(input, orcapath):
    output_file = f"{input}.out"
    with open(output_file, "w") as f:
        subprocess.run([orcapath, f"{input}.inp"], check=True, stdout=f, stderr=subprocess.STDOUT)

def process_docking_output(inp_file_path, curr_charge, curr_multiplicity, guest, host, role_key, ingredient_map, logger):
    role_name = role_key[0]
    all_results = split_docker_results(f"{inp_file_path}.docker.struc1.all.optimized.xyz", logger)

    # TODO reuse previous Result objects... eopt = new_eopt, einter = einter + new_einter
    # TODO next docking should use all poses that satisfy constraints
    c = 0
    for result in all_results:
        (path, eopt, einter, df) = result
        # TODO evaluate output for satisfying constraints
        new_path = os.path.join(os.path.dirname(path), f"struct{c}", "host.xyz")
        os.makedirs(os.path.dirname(new_path), exist_ok=True)
        shutil.move(path, new_path)
        # Edit df to assign ingredient names to atoms
        df[["ROLE", "ING", "DISH"]] = None, None, None
        n_host, n_guest = host.n_atoms, guest.n_atoms
        # transfer original labels
        for col in ["ATOM", "ATOM_ID", "ATOM_NAME", "RES_NAME", "CHAIN_ID", "RES_ID", "OCCUPANCY", "BETAFACTOR", "ELEMENT", "ROLE", "ING"]:
            df.loc[:n_host - 1, col] = host.df[col].values
            df.loc[n_host:n_host + n_guest - 1, col] = guest.df[col].values
        # transfer dish labels
        df.loc[:n_host - 1, "DISH"] = host.df["DISH"].values
        df.loc[n_host:n_host + n_guest - 1, "DISH"] = role_name

        # Make sure each molecule gets a new CHAIN_ID so that ATOM_IDs can be individually 1-based per each molecule
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
        ingredient_map[new_role_key] = product_obj
        c+=1
    return ingredient_map

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
            role_key = (str(role.title), str(guest.name), str(combo_idx)) # is a tuple
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
    prev_results = []
    prev_result_dirs = []
    for i, role in enumerate(roles):
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

            inp_file_path, curr_charge, curr_multiplicity = write_docking_input(role_key[0], role_desc.guests, role_desc.host, role_desc.constraints, workdir, qmMethod, quick, nprocs)

            # run_docking(inp_file_path, orcapath)
            inp_file_path = "/home/mchrnwsk/prometheozyme/runs/output_2025-11-03_19-19-28/base1/his/0/dock"
            ingredient_map = process_docking_output(inp_file_path, curr_charge, curr_multiplicity, role_desc.guests, role_desc.host, role_key, ingredient_map, logging)
            print(ingredient_map.items())
            [prev_results.append(ingredient_map[key]) for key in ingredient_map.keys() if key[:3] == role_key]
            print(prev_results)
            sys.exit()
 
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