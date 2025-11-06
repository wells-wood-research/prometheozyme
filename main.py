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

from define import Ingredient, Restraint, Course, col_order, col_types
from utils import split_docker_results
from evaluate import filter_conformations
from drOrca import make_orca_input
import pdbUtils

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
    # TODO are all fields of these objects useful?
    ingredients = {}
    for ing in ingredients_cfg:
        ingredient_obj = Ingredient(
            path=ing['path'],
            name=ing['name'],
            eopt=0,
            einter=0,
            charge=ing['charge'],
            multiplicity=ing['multiplicity'],
            flavours=ing.get('flavours', {}),
            df=None
        )
        ingredients[ingredient_obj.name] = ingredient_obj
    logging.debug(f"""Ingredient map is:
                {ingredients}\n""")

    # Create 'Course' objects
    # TODO are all fields of these objects useful?
    courses = {}
    for course in courses_cfg:
        if course['name'] == 'init':
            # At any given time course.host is a single value describing name of a course to take guests from for docking restraints
            # Or the one currently considered host, however keep name as plural for underlying idea
            host_candidate = None
            # There is a list of mutliple potential guests for any course
            guest_candidates = [ingredients["sub"]]
            restraints = None
        else:
            # Map candidates to Ingredient objects
            host_candidate = course['host_candidates']['name'] # this will be name of a previous Course - not actual molecular hosts yet
            guest_candidates = [ingredients[cand['name']] for cand in course['guest_candidates']]
            
            # Handle restraints if present
            restraints_data = course.get('restraints', [])
            restraints = []
            for restr in restraints_data:
                restraint = Restraint(
                    guestIdx=restr['guestIdx'],
                    hostIdx=restr['hostIdx'],
                    val=restr['val'],
                    force=restr['force']
                )
                restraints.append(restraint)
        
        # Create Course object
        course_obj = Course(
            name=course['name'],
            host=host_candidate,
            guests=guest_candidates,
            restraints=restraints
        )
        courses[course_obj.name] = course_obj 
    logging.debug(f"""Course map is:
                {courses}""")

    return courses, ingredients

def setup(configPath):
    config, outdir, orca = setup_config(configPath)
    setup_logging(outdir)
    # Prepare ingredients and courses from the configuration    
    courses, ingredients = setup_ingredients(config)
    return outdir, orca, courses, ingredients

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

########################
## DOCKING
########################

def calculate_charge(charges):
    return int(sum(charges))
    
def calculate_multiplicity(multiplicities):
    # Assume molecules are weakly or non-interacting
    spin = int(sum(m - 1 for m in multiplicities) / 2)
    multiplicity = 2*spin + 1
    return multiplicity

def dock(outdir, course_key, course_desc, orca):
    # Read orca parameters
    orcapath = orca.get("orcapth", "./orca")
    qmMethod = orca.get("qmMethod", "XTB2")
    strategy = orca.get("strategy", "normal")
    optLevel = orca.get("optLevel", "sloppyopt")
    nOpt = orca.get("nOpt", 5)
    gridExtent = orca.get("gridExtent", 15)
    nprocs = orca.get("nprocs", 8) 

    # Output of each docking step is organised as course/guest/restraints_combination
    workdir = os.path.join(outdir, *course_key)
    os.makedirs(workdir, exist_ok=True)

    inp_file_path, curr_charge, curr_multiplicity = write_docking_input(course_key[0], course_desc.guests, course_desc.host, course_desc.restraints, workdir, qmMethod, strategy, optLevel, nOpt, gridExtent, nprocs)

    run_docking(inp_file_path, orcapath)
    results_map = process_docking_output(inp_file_path, curr_charge, curr_multiplicity, course_desc.guests, course_desc.host, course_key, logging)
    return results_map

def write_docking_input(course_name, guest, host, restraints, workdir, qmMethod, strategy, optLevel, nOpt, gridExtent, nprocs):
    inp_file_path = os.path.join(workdir, f"dock")
    title = f"ORCA DOCKER: Automated Docking Algorithm for\n # course: {course_name}\n # host: {host.name}\n # guest: {guest.name}"

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
    for restraint in restraints:
        bias = {"atoms": [restraint.guestIdx, restraint.hostIdx], "val": restraint.val, "force": restraint.force}
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

def process_docking_output(inp_file_path, curr_charge, curr_multiplicity, guest, host, course_key, logger):
    course_name = course_key[0]
    all_results = split_docker_results(f"{inp_file_path}.docker.struc1.all.optimized.xyz", logger)
    results_map = {}

    # TODO reuse previous Result objects... eopt = new_eopt, einter = einter + new_einter
    # TODO next docking should use all poses that satisfy restraints
    c = 0
    for result in all_results:
        (path, eopt, einter, df) = result
        # TODO evaluate output for satisfying restraints
        # TODO after optmisation there might be duplicated results - need to remove
        new_path = os.path.join(os.path.dirname(path), f"struct{c}", "host.xyz")
        os.makedirs(os.path.dirname(new_path), exist_ok=True)
        shutil.move(path, new_path)
        # Edit df to assign ingredient names to atoms
        df[["FLAVOUR", "ING", "DISH"]] = None, None, None
        n_host, n_guest = host.n_atoms, guest.n_atoms
        # transfer original labels

        # TODO probably host df here is just the constrained element, not all
        
        for col in ["ATOM", "ATOM_ID", "ATOM_NAME", "RES_NAME", "CHAIN_ID", "RES_ID", "OCCUPANCY", "BETAFACTOR", "ELEMENT", "FLAVOUR", "ING"]:
            df.loc[:n_host - 1, col] = host.df[col].values
            df.loc[n_host:n_host + n_guest - 1, col] = guest.df[col].values
        # transfer dish labels
        df.loc[:n_host - 1, "DISH"] = host.df["DISH"].values
        df.loc[n_host:n_host + n_guest - 1, "DISH"] = course_name

        # Make sure each molecule gets a new CHAIN_ID so that ATOM_IDs can be individually 1-based per each molecule
        print(df)
        df = assign_chain_ids(df)

        # TODO ATOM_ID should probably continue... of chain id should be changed...
        df = df[col_order].astype(col_types)

        # Save PDB
        df_to_save = df.loc[:, ~df.columns.isin(["FLAVOUR", "ING", "DISH"])]
        pdbUtils.df2pdb(df_to_save, new_path.replace(".xyz", ".pdb"))
        print(f"Saved to {new_path.replace('.xyz', '.pdb')}")
        
        new_course_key = course_key + (str(c), )
        product_obj = Ingredient(
            path=new_path,
            name=("_").join(new_course_key),
            eopt=eopt,
            einter=einter,
            charge=curr_charge,
            multiplicity=curr_multiplicity,
            df=df
        )
        results_map[new_course_key] = product_obj
        c+=1
    return results_map

########################
## RESTRAINTS
########################

def assign_restraint_idx(guest, host, restraint):
    # TODO add error prints here
    guest_matches = [(i, row["ATOM_NAME"], row["FLAVOUR"]) for i, row in guest.df.iterrows() if restraint.guestIdx in row["FLAVOUR"]]
    logging.debug("Expanding restraints for guest:")
    logging.debug(guest.df)
    logging.debug(f"Original guestIdx to be restrained: {restraint.guestIdx}")
    logging.debug(f"Actual guest atoms matching the original: {guest_matches}\n")
    
    host_matches = [(i, row["ATOM_NAME"], row["FLAVOUR"]) for i, row in host.df.iterrows() if restraint.hostIdx in row["FLAVOUR"]]
    logging.debug("Expanding restraints for host:")
    logging.debug(host.df)
    logging.debug(f"Original hostIdx to be restrained: {restraint.hostIdx}")
    logging.debug(f"Actual host atoms matching the original: {host_matches}\n")

    return guest_matches, host_matches

def expand_restraints(guest, host, restraint):
    # Filter tuples whose activity matches the requested flavour name
    guest_matches, host_matches = assign_restraint_idx(guest, host, restraint)

    # To define bond bias potential in ORCA DOCKING step
    bias_params = []
    for g_item, h_item in itertools.product(guest_matches, host_matches):
        # g_item and h_item are tuples like (index, name, activity)
        g_idx = g_item[0]  # numeric index in the PDB/molecule
        h_idx = h_item[0]
        bias_params.append((g_idx, h_idx, restraint.val, restraint.force))
    return bias_params

def expand_ingredient_and_restraint_combinations(course):
    expanded_course = {}

    # Loop all candidate host/guest pairs
    host = course.host
    for guest in course.guests:
        # For each restraint compute the list of concrete options (g_idx, h_idx, val, force)
        bias_params_per_restraint = []
        # TODO how would this work for multiple course restraints, e.g. distance and angle?
        # TODO desired behaviour is to make all combinations of one restraint, and all combinations of the other (separately), and then combine them directlyin all possibilities
        for restraint in course.restraints or []:
            bias_params = expand_restraints(guest, host, restraint)
            # If a restraint has no matching atoms, this course/host/guest pair can't satisfy it:
            if not bias_params:
                # skip this host/guest pair entirely
                bias_params_per_restraint = []
                break
            bias_params_per_restraint.append(bias_params)

        # Case 1) no restraints were defined for this course:
        if not course.restraints:
            # Create single course copy with no concrete restraints
            new_course = copy.deepcopy(course)
            new_course.host = host
            new_course.guests = guest
            new_course.restraints = []
            course_key = f"{course.name}_{host.name}_{guest.name}_constrX"
            logging.debug(f"No restraints defined for course {course_key}")
            expanded_course[course_key] = new_course
            continue

        # Case 2) restraints cannot be satisfied
        if not bias_params_per_restraint:
            # TODO add error prints here
            continue

        # Cartesian product across lists of options (one option per restraint)
        # If one ingredient has multiple atoms of same flavour (e.g. multiple h_donor),
        # this creates multiple combinations (a "mash-up"), one for each of the atoms
        for mash_idx, mash in enumerate(itertools.product(*bias_params_per_restraint)):
            # 'mash' is a tuple of option tuples, one per restraint
            # Build new restraint tuples representing these concrete restraint
            mashed_restraints = []
            for orig_restraint, bias_params in zip(course.restraints, mash):
                g_idx, h_idx, val, force = bias_params
                # Create a copy of the original restraint but with resolved atom indices
                new_restraint = copy.deepcopy(orig_restraint)
                # Overwrite guestIdx/hostIdx with the concrete atom indices
                new_restraint.guestIdx = g_idx
                new_restraint.hostIdx = h_idx
                # Ensure val/force are set to the chosen option (val likely same as orig)
                new_restraint.val = val
                new_restraint.force = force
                mashed_restraints.append(new_restraint)

            # Make new course copy and set host/guest to the concrete candidates and restraints
            new_course = copy.deepcopy(course)
            new_course.host = host
            new_course.guests = guest  # TODO guests vs guest... typing... should this be a list here?
            new_course.restraints = mashed_restraints

            # Unique key: include the combination index so expansion products of same host and guest don't overwrite themselves
            # is a tuple
            course_key = (str(course.name), str(guest.name), str(mash_idx))
            expanded_course[course_key] = new_course

    return expanded_course

########################
## MAIN LOOP
########################

def main(args):
    # Read config file
    outdir, orca, courses, ingredients = setup(args.config)

    # Results of previous course, carried over as starting point (HOST) of next docking step
    leftovers = []
    
    # Prepare the substrate as a fake result of first docking for processing consistency
    init_host = courses["init"].guests[0]
    init_host_dir = os.path.join(outdir, f"base_stock")
    init_host_path = os.path.join(init_host_dir, "host.pdb") # TODO do you want "leftovers" to be PDB or XYZ based? DO you want to always (!) write both XYZ and PDB?
    os.makedirs(init_host_dir)
    shutil.copy(init_host.path, init_host_path)
    # TODO assumes there's xyz present - need to expect user only to prepare PDBs
    shutil.copy(init_host.path.replace(".pdb", ".xyz"), init_host_path.replace(".pdb", ".xyz"))
    init_host.path = init_host_path
    leftovers.append(init_host)

    # Run ORCA DOCKER for each course, exhaustive for host/guest mix
    for i, (course_name, course) in enumerate(courses.items()):
        if course_name == "init":
            continue
        for serving in leftovers:
            # There's two kinds of "HOSTS" in every docking step:
            # 1) host file which is the complete product of previous course
            # 2) host molecule to restrain guest towards (for seleciton of atom indices in BIAS block of ORCA input file)
            new_host_for_docking_ing = serving
            new_host_for_docking_df = new_host_for_docking_ing.df
            new_host_for_restraints_df = new_host_for_docking_df[new_host_for_docking_df["DISH"] == course.host]
            new_host_for_restraints_name = new_host_for_restraints_df["ING"].unique()[0]
            new_host_for_restraints_ing = ingredients[new_host_for_restraints_name]
            new_host_for_restraints_ing.df = new_host_for_restraints_df

            new_course_for_restraints = copy.deepcopy(course)
            new_course_for_restraints.host = new_host_for_restraints_ing

            # TODO I am passing here a fake host for restraints which e..g is only df of molecule C from ABCD
            # TODO does that interrupt with other things?
            expanded_course = expand_ingredient_and_restraint_combinations(new_course_for_restraints)
            logging.debug(f"""Expanded ingredient and restraint combinations for
                        course name: {course.name}
                        number of combinations: {len(expanded_course)}
                        combination keys: {list(expanded_course.keys())}""")

            new_leftovers = []
            for course_key, course_desc in expanded_course.items():
                logging.info(f"Processing course: {course_key}")
                course_desc.host = new_host_for_docking_ing
                # TODO don't nest directories so deep - more vertical structure
                outdir = os.path.dirname(serving.path)
                waste_bucket = dock(outdir, course_key, course_desc, orca)
                [new_leftovers.append(waste_bucket[key]) for key in waste_bucket.keys() if key[:3] == course_key]
        leftovers = new_leftovers

    logging.info("""
          Main script complete!
          ⚝⭒٭⋆⚝⭒٭⋆⚝⭒٭⋆⚝⭒٭⋆⚝⭒٭⋆⚝⭒٭⋆
          """)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Process some ingredients and courses.")
    parser.add_argument('--config', type=str, default=None, help="Path to the configuration YAML file")
    args = parser.parse_args()
    if not args.config:
        # TODO del static assignemnt when finished testing
        args.config = "/home/mchrnwsk/prometheozyme/config.yaml"
    main(args)