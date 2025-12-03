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
import uuid

from utils.drThing import Ingredient, Restraint, Course, col_order, col_types
from utils.drStructure import extract_ok_docker_results, isPDB, isXYZ, pdb2df, df2pdb, df2xyz, pdb2xyz, xyz2df, xyz2pdb, write_multi_pdb
from utils.drOrca import make_orca_input

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
    verbosity = misc.get("verbosity", ".")

    # Get orca docker parameters
    orca = config.get("orca", {})

    # Setup output dir
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    outdir = os.path.join(workdir, f"output_{timestamp}")
    os.makedirs(outdir, exist_ok=True)

    return config, outdir, orca, verbosity

def setup_logging(workdir, verbosity="info"):
    # Define a custom logging level called VERBOSE
    VERBOSE_LEVEL_NUM = 5
    logging.addLevelName(VERBOSE_LEVEL_NUM, "VERBOSE")
    def verbose_root(message, *args, **kwargs):
        logging.log(VERBOSE_LEVEL_NUM, message, *args, **kwargs)
    logging.verbose = verbose_root

    # Map string verbosity to logging levels
    level_map = {
        "verbose": VERBOSE_LEVEL_NUM,
        "debug": logging.DEBUG,
        "info": logging.INFO,
        "warning": logging.WARNING,
        "error": logging.ERROR,
        "critical": logging.CRITICAL,
    }
    log_level = level_map.get(verbosity.lower(), logging.INFO)
    
    # Remove existing handlers
    root_logger = logging.getLogger()
    for handler in list(root_logger.handlers):
        root_logger.removeHandler(handler)

    # Log to both console and file
    handlers = [logging.StreamHandler()]  # console
    handlers.append(logging.FileHandler(os.path.join(workdir, "log.txt")))  # file

    logging.basicConfig(
        format="%(asctime)s - %(levelname)s - %(message)s",
        level=log_level,
        handlers=handlers,
    )

def setup_ingredients(config):
    allOk = True
    ingredients_cfg = config.get("ingredients", [])
    courses_cfg = config.get("courses", [])

    logging.debug(f"""Ingredients are:
                {ingredients_cfg}\n""")
    logging.debug(f"""Courses are:
                {courses_cfg}\n""")

    # Create 'Ingredient' objects
    ingredients = {}
    for ing in ingredients_cfg:
        pathPDB=ing['path']
        pathXYZ=pathPDB.replace(".pdb", ".xyz")
        pdb2xyz(pathPDB=pathPDB, outXYZ=pathXYZ, logger=logging)
        ingredient_obj = Ingredient(
            pathPDB=pathPDB,
            pathXYZ=pathXYZ,
            name=ing['name'],
            eopt=0,
            einter=0,
            charge=ing['charge'],
            multiplicity=ing['multiplicity'],
            flavours=ing.get('flavours', {}),
            df=None
        )
        ingredients[ingredient_obj.name] = ingredient_obj
    # TODO to print this properly, add some "print" function to the classes
    # (otherwise prints things like: {'sub': <utils.drThing.Ingredient object at 0x77104e150290> ... or ... Restraints are: [<utils.drThing.Restraint object at 0x77104ea16490>])
    logging.debug(f"""Ingredient map is:
                {ingredients}\n""")

    # Create 'Course' objects
    courses = {}
    for course in courses_cfg:
        if course['name'] == 'init':
            # At any given time course.host is a single value describing name of a course to take guests from for docking restraints
            # Or the one currently considered host, however keep name as plural for underlying idea
            host_candidate = None
            # There is a list of mutliple potential guests for any course
            guest_candidates = [ingredients[cand['name']] for cand in course['guest_candidates']]
            restraints = None
            orcaSettings = None
        else:
            # Map candidates to Ingredient objects
            host_candidate = course['host_candidates']['name'] # this will be name of a previous Course - not actual molecular hosts yet
            guest_candidates = [ingredients[cand['name']] for cand in course['guest_candidates']]
            
            # Handle restraints if present
            restraints_data = course.get('restraints', [])
            restraints = []
            for restr in restraints_data:
                restraint = Restraint(
                    guestIdx=restr.get('guestIdx', None),
                    hostIdx=restr.get('hostIdx', None),
                    val=restr.get('val', None),
                    tol=restr.get('tol', 0.5),
                    force=restr.get('force', 100)
                )
                if restraint.guestIdx is None or restraint.hostIdx is None or restraint.val is None:
                    logging.warning(f"No guestIdx, hostIdx, val defined for course {course['name']} restraints. Fix or remove the restraints block from config file.")
                    allOk = False
                    restraint = None 
                restraints.append(restraint)
        
        # Create Course object
        course_obj = Course(
            name=course.get("name", str(uuid.uuid4())),
            host=host_candidate,
            guests=guest_candidates,
            restraints=restraints,
            orcaSettings=course.get("orcaSettings", None)
        )
        courses[course_obj.name] = course_obj 
        logging.debug(f"Restraints are: {course_obj.restraints}")
    logging.debug(f"""Course map is:
                {courses}\n""")

    return courses, ingredients, allOk

def setup(configPath):
    allOk = True
    config, outdir, orca, verbosity = setup_config(configPath)
    setup_logging(outdir, verbosity)
    # Prepare ingredients and courses from the configuration    
    courses, ingredients, is_step_ok = setup_ingredients(config)
    if not is_step_ok:
        allOk = False
    return outdir, orca, courses, ingredients, allOk

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

def define_orca_params(orca, course_desc):
    # Define default ORCA settings
    default_orca_settings = {
        "orcapth": "./orca",
        "qmMethod": "XTB2",
        "strategy": "NORMAL",
        "optLevel": "sloppyopt",
        "nOpt": 5,
        "fixHost": True,
        "gridExtent": 15,
        "nprocs": 8
    }

    # Merge global ORCA config (orca dict from config file)
    # This overrides defaults with general config-level values
    merged_orca_settings = {**default_orca_settings, **orca}

    # Check if this course defines specific ORCA overrides
    course_orca_settings = course_desc.orcaSettings
    if course_orca_settings:
        # Only overwrite recognized keys
        for key in default_orca_settings.keys():
            if key in course_orca_settings:
                merged_orca_settings[key] = course_orca_settings[key]

    # Extract final ORCA parameters
    orcapath   = merged_orca_settings["orcapth"]
    qmMethod   = merged_orca_settings["qmMethod"]
    strategy   = merged_orca_settings["strategy"]
    optLevel   = merged_orca_settings["optLevel"]
    nOpt       = merged_orca_settings["nOpt"]
    fixHost    = merged_orca_settings["fixHost"]
    gridExtent = merged_orca_settings["gridExtent"]
    nprocs     = merged_orca_settings["nprocs"]

    return orcapath, qmMethod, strategy, optLevel, nOpt, fixHost, gridExtent, nprocs

def dock(outdir, course_key, course_desc, orca):
    try:
        orcapath, qmMethod, strategy, optLevel, nOpt, fixHost, gridExtent, nprocs = define_orca_params(orca, course_desc)
        
        workdir = os.path.join(outdir, *course_key[1:])
        os.makedirs(workdir, exist_ok=True)

        biases = process_biases(course_desc.restraints)
        inp_file_path, curr_charge, curr_multiplicity = write_docking_input(
            course_key[0], course_desc.guests[0], course_desc.host, biases,
            workdir, qmMethod, strategy, optLevel, nOpt, fixHost, gridExtent, nprocs
        )
        logging.info(f"Docking input written at path: {inp_file_path}.inp")
        logging.info(f"Running docking...")
        run_docking(inp_file_path, orcapath)
        logging.info(f"Docking complete. See details at path: {inp_file_path}.out")

        results_map = process_docking_output(
            inp_file_path, curr_charge, curr_multiplicity,
            course_desc.guests[0], course_desc.host, biases, course_key
        )
        if not results_map:
            logging.warning(f"Docking at {os.path.dirname(inp_file_path)} produced no usable results.")

        return results_map

    except subprocess.CalledProcessError as e:
        logging.exception(f"Critical error during ORCA run: {e}")
        return {}
    except Exception as e:
        logging.exception(f"Unexpected error in docking step: {e}")
        return {}

def process_biases(restraints):
    # Define bond bias potential to impose restraints of distance between guest and host atoms
    # Must be defined by atom numbers for GUEST FIRST, then HOST SECOND, absolute indexing for each molecule independently
    # For example, to add a bond bias between atom 2 from the GUEST and atom 19 from the HOST, write: BIAS { B 2 19 } END 
    biases = []
    for restraint in restraints:
        bias = {"atoms": [restraint.guestIdx, restraint.hostIdx], "val": restraint.val, "tol": restraint.tol, "force": restraint.force}
        biases.append(bias)
    return biases

def write_docking_input(course_name, guest, host, biases, workdir, qmMethod, strategy, optLevel, nOpt, fixHost, gridExtent, nprocs):
    inp_file_path = os.path.join(workdir, f"dock")
    title = f"ORCA DOCKER: Automated Docking Algorithm for\n # course: {course_name}\n # host: {host.name}\n # guest: {guest.name}"

    # charge and multiplicity of host only - guest is defined under %DOCKER GUESTCHARGE and GUESTMULT
    moleculeInfo = {"charge": host.charge, "multiplicity": host.multiplicity}

    docker = {"guestPath": guest.pathXYZ, "guestCharge": guest.charge, "guestMultiplicity": guest.multiplicity, "fixHost": fixHost, "bias": biases, "strategy": strategy, "optLevel": optLevel, "nOpt": nOpt, "gridExtent": gridExtent}

    # Use the updated XYZ file for optimization
    make_orca_input(orcaInput=inp_file_path,
                    title=title,
                    simpleInputLine=[qmMethod],
                    inputFormat="xyzfile",
                    inputFile=host.pathXYZ,
                    moleculeInfo=moleculeInfo,
                    parallelize=nprocs,
                    docker=docker)
    
    # Metadata returned to facilitate further docking, reusing products of earlier docking as hosts
    return inp_file_path, calculate_charge([guest.charge, host.charge]), calculate_multiplicity([guest.multiplicity, host.multiplicity])

def run_docking(input, orcapath):
    stdout_file = f"{input}.out"
    stderr_file = f"{input}.err"
    with open(stdout_file, "w") as out, open(stderr_file, "w") as err:
        result = subprocess.run([orcapath, f"{input}.inp"], check=True, stdout=out, stderr=err)
    
    if result.returncode != 0:
        with stdout_file.open("r", errors="ignore") as f:
            lines = f.readlines()
        for l in lines[-13:]:
            print(l.rstrip())

def process_docking_output(inp_file_path, curr_charge, curr_multiplicity, guest, host, biases, course_key):
    course_name = course_key[0]
    # The following function splits mutli-xyz output of docker,
    # evaluates each result for satisfaction of constraints,
    # and writes only the correct outputs to separate single-xyz files
    ok_results = extract_ok_docker_results(f"{inp_file_path}.docker.struc1.all.optimized.xyz", host.n_atoms, biases, logger=logging)
    results_map = {}

    c = 0
    for result in ok_results:
        (pathXYZ, eopt, einter, df) = result
        new_pathXYZ = os.path.join(os.path.dirname(pathXYZ), f"serving{c}", "host.xyz")
        os.makedirs(os.path.dirname(new_pathXYZ), exist_ok=True)
        shutil.move(pathXYZ, new_pathXYZ)
        # Edit df to assign ingredient names to atoms
        df[["FLAVOUR", "ING", "DISH"]] = None, None, None
        n_host, n_guest = host.n_atoms, guest.n_atoms

        for col in ["ATOM", "ATOM_ID", "ATOM_NAME", "RES_NAME", "CHAIN_ID", "RES_ID", "OCCUPANCY", "BETAFACTOR", "ELEMENT", "FLAVOUR", "ING"]:
            df.loc[:n_host - 1, col] = host.df[col].values
            df.loc[n_host:n_host + n_guest - 1, col] = guest.df[col].values
        # transfer dish labels
        df.loc[:n_host - 1, "DISH"] = host.df["DISH"].values
        df.loc[n_host:n_host + n_guest - 1, "DISH"] = course_name

        # Make sure each molecule gets a new CHAIN_ID so that ATOM_IDs can be individually 1-based per each molecule
        df = assign_chain_ids(df)
        df = df[col_order].astype(col_types)
        logging.verbose( f"Result of this docking step:\n{df}")

        # Save PDB
        df_to_save = df.loc[:, ~df.columns.isin(["FLAVOUR", "ING", "DISH"])]
        new_pathPDB = new_pathXYZ.replace(".xyz", ".pdb")
        df2pdb(df=df_to_save, outPDB=new_pathPDB, logger=logging)
        logging.info(f"Saved result to {new_pathPDB}")
        
        new_course_key = course_key + (str(c), )
        product_obj = Ingredient(
            pathPDB=new_pathPDB,
            pathXYZ=new_pathXYZ,
            name=f"result of {new_course_key}",
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

def assign_restraint_idx(guest, host, restraint, course_name):
    guest_matches = [(i, row["ATOM_NAME"], row["FLAVOUR"]) for i, row in guest.df.iterrows() if restraint.guestIdx in row["FLAVOUR"]]
    if not guest_matches:
        logging.error(f"Guest {guest.name} has no atoms matching flavour of restraint {restraint.guestIdx} for course {course_name}\n")
    else:
        logging.verbose("Expanding restraints for guest:")
        logging.verbose(f"\n{guest.df}")
        logging.debug(f"Original guestIdx to be restrained: {restraint.guestIdx}")
        logging.debug(f"Guest atoms matching the restraint (idx, ATOM_NAME, FLAVOUR): {guest_matches}")
    
    host_matches = [(i, row["ATOM_NAME"], row["FLAVOUR"]) for i, row in host.df.iterrows() if restraint.hostIdx in row["FLAVOUR"]]
    if not host_matches:
        logging.error(f"Host {host.name} has no atoms matching flavour of restraint {restraint.hostIdx} for course {course_name}")
    else:
        logging.verbose("Expanding restraints for host:")
        logging.verbose(f"\n{host.df}")
        logging.debug(f"Original hostIdx to be restrained: {restraint.hostIdx}")
        logging.debug(f"Host atoms matching the restraint (idx, ATOM_NAME, FLAVOUR): {host_matches}")

    return guest_matches, host_matches

def expand_restraints(guest, host, restraint, course_name):
    # Filter tuples whose activity matches the requested flavour name
    guest_matches, host_matches = assign_restraint_idx(guest, host, restraint, course_name)

    # To define bond bias potential in ORCA DOCKING step
    bias_params = []
    for g_item, h_item in itertools.product(guest_matches, host_matches):
        # g_item and h_item are tuples like (index, name, activity)
        g_idx = g_item[0]  # numeric index in the PDB/molecule
        h_idx = h_item[0]
        bias_params.append((g_idx, h_idx, restraint.val, restraint.tol, restraint.force))
    return bias_params

def expand_ingredient_and_restraint_combinations(course, host):
    allOk = True
    expanded_course = {}
    # The desired behaviour is to make all combinations of one restraint, and all combinations of the other (separately), and then combine them directly in all possibilities
    # The code loops over all candidate guest-host ingredient pairs,
    # for each restraint, expands to all concrete atom index pairs,
    # for multiple restraints, takes full Cartesian product,
    # and produces a new Course object per combination

    # Loop all candidate host/guest pairs
    for guest in course.guests:
        # For each restraint compute the list of concrete options (g_idx, h_idx, val, tol, force)
        bias_params_per_restraint = []
        for restraint in course.restraints or []:
            # from setup_ingredients, course.restraints could be an empty list or a list on Nones
            if restraint is None:
                continue
            bias_params = expand_restraints(guest, host, restraint, course.name)
            # If a restraint has no matching atoms, this course/host/guest pair can't satisfy it:
            if not bias_params:
                logging.error(f"Course {course.name} has restraints defined (e.g. guestIdx {course.restraints[0].guestIdx}, hostIdx {course.restraints[0].hostIdx}, val: {course.restraints[0].val}) but bias potentials could not be defined.")
                allOk = False
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
            logging.info(f"No restraints defined for course {course_key}\n")
            expanded_course[course_key] = new_course
            continue

        # Case 2) restraints cannot be satisfied
        if not bias_params_per_restraint:
            logging.error(f"No bias params could be defined for restraints of course {course.name}.")
            allOk = False
            break

        # Cartesian product across lists of options (one option per restraint)
        # If one ingredient has multiple atoms of same flavour (e.g. multiple h_donor),
        # this creates multiple combinations (a "mash-up"), one for each of the atoms
        for mash_idx, mash in enumerate(itertools.product(*bias_params_per_restraint)):
            # 'mash' is a tuple of option tuples, one per restraint
            # Build new restraint tuples representing these concrete restraint
            mashed_restraints = []
            for orig_restraint, bias_params in zip(course.restraints, mash):
                g_idx, h_idx, val, tol, force = bias_params
                # Create a copy of the original restraint but with resolved atom indices
                new_restraint = copy.deepcopy(orig_restraint)
                # Overwrite guestIdx/hostIdx with the concrete atom indices
                new_restraint.guestIdx = g_idx
                new_restraint.hostIdx = h_idx
                # Ensure val/tol/force are set to the chosen option
                new_restraint.val = val
                new_restraint.tol = tol
                new_restraint.force = force
                mashed_restraints.append(new_restraint)

            # Make new course copy and set host/guest to the concrete candidates and restraints
            new_course = copy.deepcopy(course)
            new_course.host = host
            new_course.guests = [guest]
            new_course.restraints = mashed_restraints

            # Unique key: include the combination index so expansion products of same host and guest don't overwrite themselves
            course_key = (str(course.name), str(guest.name), str(mash_idx))
            expanded_course[course_key] = new_course

    return expanded_course, allOk

########################
## MAIN LOOP
########################

def main(args):
    allOk = True
    # Read config file
    outdir, orca, courses, ingredients, allOk = setup(args.config)

    # Copy config file to output directory for reproducibility / tracking
    shutil.copy(args.config, outdir)
    
    logging.info(f"Cooking begins - recipe from {args.config}")

    # Results of previous course, carried over as starting point (HOST) of next docking step
    leftovers = []
    
    # Prepare the substrate as a fake result of first docking for processing consistency
    try:
        init_host = courses["init"].guests[0]
        init_host_dir = os.path.join(outdir, f"base_stock")
        init_host_pathPDB = os.path.join(init_host_dir, "host.pdb")
        init_host_pathXYZ = init_host_pathPDB.replace(".pdb", ".xyz")
        os.makedirs(init_host_dir)
        shutil.copy(init_host.pathPDB, init_host_pathPDB)
        shutil.copy(init_host.pathXYZ, init_host_pathXYZ)
        init_host.pathPDB = init_host_pathPDB
        init_host.pathXYZ = init_host_pathXYZ
        leftovers.append(init_host)
    except Exception as e:
        logging.exception(f"Critical error while preparing base stock (initial host) from substrate: {e}")
        allOk = False
        return

    # Run ORCA DOCKER for each course, exhaustive for host/guest mix
    for i, (course_name, course) in enumerate(courses.items()):
        if not allOk:
            logging.error("Terminating remaining docking steps due to previous critical error.")
            break

        if course_name == "init":
            continue

        new_leftovers = []
        for j, serving in enumerate(leftovers):
            try:
                # There's two kinds of "HOSTS" in every docking step:
                # 1) host file which is the complete product of previous course
                # 2) host molecule to restrain guest towards (for seleciton of atom indices in BIAS block of ORCA input file)
                new_host_for_docking_ing = serving
                new_host_for_docking_df = new_host_for_docking_ing.df
                new_host_for_restraints_df = new_host_for_docking_df[new_host_for_docking_df["DISH"] == course.host]
                new_host_for_restraints_name = new_host_for_restraints_df["ING"].unique()[0]
                new_host_for_restraints_ing = ingredients[new_host_for_restraints_name]
                new_host_for_restraints_ing.df = new_host_for_restraints_df
            except Exception as e:
                logging.exception(f"Critical error preparing host for course {course_name}: {e}")
                allOk = False
                break

            expanded_course, allOk = expand_ingredient_and_restraint_combinations(course, new_host_for_restraints_ing)
            logging.info(f"""Expanded ingredient and restraint combinations for course {course.name} ({len(expanded_course)} combinations),
                            keys: {list(expanded_course.keys())}""")
            
            for course_key, course_desc in expanded_course.items():
                try:
                    logging.info(f"""
                                Processing course: {course_key}""")
                    course_desc.host = new_host_for_docking_ing
                    curr_outdir = os.path.join(outdir, f"course{i}_{course_key[0]}", f"serving{j}")
                    waste_bucket = dock(curr_outdir, course_key, course_desc, orca)
                    [new_leftovers.append(waste_bucket[key]) for key in waste_bucket.keys() if key[:3] == course_key]
                except Exception as e:
                    logging.exception(f"Critical error in course {course_key}: {e}")
                    allOk = False
                    break
        leftovers = new_leftovers
        if len(leftovers) < 1:
            allOk = False
            logging.warning(f"Critical failure during docking of course {course_key} - docking produced no usable results.")
            break

    # Report here final theozymes
    # Sort by einter value (ascending)
    leftovers_sorted = sorted(leftovers, key=lambda x: x.einter)

    # Create output directory
    results_dir = os.path.join(outdir, "results")
    os.makedirs(results_dir, exist_ok=True)

    logging.info("Successfully cooked the following theozymes:")

    # Copy files with renamed names based on order
    for i, ing in enumerate(leftovers_sorted, start=1):
        new_name = f"result{i}"
        new_pathPDB = os.path.join(results_dir, f"{new_name}.pdb")
        # Copy and rename files
        shutil.copy(ing.pathPDB, new_pathPDB)
        # Log the results
        logging.info(
            f"path: {new_pathPDB}, eopt: {ing.eopt} (Eh), einter: {ing.einter} (kcal/mol)"
        )
    
    # Merge into one multi-frame PDB file for easier analysis
    result_paths = [os.path.join(results_dir, x) for x in os.listdir(results_dir) if x.lower().endswith(".pdb")]
    # TODO merging duplicates atoms
    write_multi_pdb(result_paths, os.path.join(results_dir, "merged.pdb"))

    if allOk:
        logging.info("""
            Cooking complete - bon appetit!
            ⚝⭒٭⋆⚝⭒٭⋆⚝⭒٭⋆⚝⭒٭⋆⚝⭒٭⋆⚝⭒٭⋆
            """)
    else:
        logging.info("""
            Oh, crepe - something's gone wrong.
                     Try different settings?
            ☠︎︎⭒٭⋆☠︎︎⭒٭⋆☠︎︎⭒٭⋆☠︎︎⭒٭⋆☠︎︎⭒٭⋆☠︎︎⭒٭⋆
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