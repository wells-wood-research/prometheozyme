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
from pathlib import Path
import rmsd
from collections import defaultdict
leftovers_dict = defaultdict(list)

from pathlib import Path
cwd = Path(__file__).resolve().parent

from utils import types, parse, workflow, structure, orca, validate, lookup

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

def populate_ingredient_metadata(ingredients, metadata_map):
    allOk = True

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
        course_name = course.get("name", str(uuid.uuid4()))
        if course_name == 'init':
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
                restraint_property = restr.get('property', 'distance')
                restraint_selection = restr.get('selection', [])
                restraint_params = restr.get('parameters', {"val": None, "uptol": None, "downtol": None, "force": None})

                restraint = Restraint(
                    property=restraint_property,
                    sele=restraint_selection,
                    params=restraint_params,
                    step=course_name
                )
                if restraint.property is None or restraint.sele is None or restraint.params is None:
                    logging.warning(f"No property, atom selection or restraint parameters defined for course {course['name']} restraints. Fix or remove the restraints block from config file.")
                    allOk = False
                    restraint = None 
                restraints.append(restraint)
        
        # Create Course object
        course_obj = Course(
            name=course_name,
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
    config, outdir, orca, verbosity, rmsd_threshold = setup_config(configPath)
    setup_logging(outdir, verbosity)
    # Prepare ingredients and courses from the configuration    
    courses, ingredients, is_step_ok = setup_ingredients(config)
    if not is_step_ok:
        allOk = False
    return outdir, rmsd_threshold, orca, courses, ingredients, allOk

# ----------------------------
# HELPER FUNCTIONS
# ----------------------------

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

# ----------------------------
# DOCKING
# ----------------------------

def calculate_charge(charges):
    return int(sum(charges))
    
def calculate_multiplicity(multiplicities):
    # Assume molecules are weakly or non-interacting
    spin = int(sum(m - 1 for m in multiplicities) / 2)
    multiplicity = 2*spin + 1
    return multiplicity

def define_orca_params(orca, course_desc=None):
    # Define default ORCA settings
    default_orca_settings = {
        "orcapath": "./orca",
        "qmMethod_dock": "XTB2",
        "qmMethod_opt": "XTB2",
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
    if course_desc:
        course_orca_settings = course_desc.orcaSettings
        if course_orca_settings:
            # Only overwrite recognized keys
            for key in default_orca_settings.keys():
                if key in course_orca_settings:
                    merged_orca_settings[key] = course_orca_settings[key]

    # Extract final ORCA parameters
    orcapath = merged_orca_settings["orcapath"]
    qmMethod_dock = merged_orca_settings["qmMethod_dock"]
    qmMethod_opt = merged_orca_settings["qmMethod_opt"]
    strategy = merged_orca_settings["strategy"]
    optLevel = merged_orca_settings["optLevel"]
    nOpt = merged_orca_settings["nOpt"]
    fixHost = merged_orca_settings["fixHost"]
    gridExtent = merged_orca_settings["gridExtent"]
    nprocs = merged_orca_settings["nprocs"]

    return orcapath, qmMethod_dock, qmMethod_opt, strategy, optLevel, nOpt, fixHost, gridExtent, nprocs

def dock(host_input_file, guest_input_file, new_restraints, old_restraints, output_dir, orca):
    return None

def dock_old(outdir, course_key, course_desc, orca):
    try:
        orcapath, qmMethod_dock, qmMethod_opt, strategy, optLevel, nOpt, fixHost, gridExtent, nprocs = define_orca_params(orca, course_desc)
        
        workdir = os.path.join(outdir, *course_key[1:])
        os.makedirs(workdir, exist_ok=True)

        host_n_atoms = course_desc.host.n_atoms
        current_step_restraints = process_restraints_for_docking(course_desc.restraints, host_n_atoms)
        restraints_abs = course_desc.host.restraints + current_step_restraints

        inp_file_path, curr_charge, curr_multiplicity = write_docking_input(
            course_key[0], course_desc.guests[0], course_desc.host, current_step_restraints,
            workdir, qmMethod_dock, strategy, optLevel, nOpt, fixHost, gridExtent, nprocs
        )
        logging.info(f"Docking input written at path: {inp_file_path}.inp")
        logging.info(f"Running docking...")
        result = run_orca(inp_file_path, orcapath, timeout=None)
        logging.info(f"Docking complete. See details at path: {inp_file_path}.out")
        if result != 0 :
            logging.error(f"Docking returned status code {result}. Skip to next theozyme...")
            return {} # TODO is this the right error handling?


        results_map = process_docking_output(
            inp_file_path, curr_charge, curr_multiplicity,
            course_desc.guests[0], course_desc.host, restraints_abs, course_key
        )
        if not results_map:
            logging.warning(f"Docking at {os.path.dirname(inp_file_path)} produced no usable results.")
            return {}
        
        return results_map

    except subprocess.CalledProcessError as e:
        logging.exception(f"Critical error during ORCA run: {e}")
        return {}
    except Exception as e:
        logging.exception(f"Unexpected error in docking step: {e}")
        return {}

def process_restraints_for_docking(restraints, host_n_atoms):
    """
    Input:  list of Restraint objects with parent='host'/'guest', idx still guest-local.
    Output: list of Restraint objects with absolute indices.
    """
    restraints_abs = []
    for r in restraints:
        r_new = copy.deepcopy(r)
        r_new.sele = [
            Selection(s.parent, _abs_index(s, host_n_atoms))
            for s in r.sele
        ]
        restraints_abs.append(r_new)
    return restraints_abs

def write_docking_input(course_name, guest, host, restraints_abs, workdir, qmMethod, strategy, optLevel, nOpt, fixHost, gridExtent, nprocs):
    inp_file_path = os.path.join(workdir, f"dock")
    title = f"ORCA DOCKER: Automated Docking Algorithm for\n # course: {course_name}\n # host: {host.name}\n # guest: {guest.name}"

    # charge and multiplicity of host only - guest is defined under %DOCKER GUESTCHARGE and GUESTMULT
    moleculeInfo = {"charge": host.charge, "multiplicity": host.multiplicity}
    
    docker_biases = []
    geom_biases = []
    for r in restraints_abs:
        if r.property == "distance":
            a_abs = r.sele[0].idx
            b_abs = r.sele[1].idx
            p = [s.parent for s in r.sele]

            if set(p) == {"guest", "host"}:
                # ensure guest first

                if r.sele[0].parent == "host":
                    a_abs, b_abs = b_abs, a_abs
                a_rel = a_abs - host.n_atoms
                docker_biases.append({
                    "atoms": [a_rel,b_abs],
                    "val": r.params.val,
                    "uptol": r.params.uptol,
                    "downtol": r.params.downtol,
                    "force": r.params.force
                })
            elif set(p) == {"host"}:
                geom_biases.append({
                    "atoms": [a_abs,b_abs],
                    "val": r.params.val,
                    "uptol": r.params.uptol,
                    "downtol": r.params.downtol,
                })
    
    docker = {"guestPath": guest.pathXYZ, "guestCharge": guest.charge, "guestMultiplicity": guest.multiplicity, "fixHost": fixHost, "bias": docker_biases, "strategy": strategy, "optLevel": optLevel, "nOpt": nOpt, "gridExtent": gridExtent}
    geom = None if not geom_biases else {"keep": geom_biases}

    # Use the updated XYZ file for optimization
    make_orca_input(orcaInput=inp_file_path,
                    title=title,
                    simpleInputLine=[qmMethod],
                    inputFormat="xyzfile",
                    inputFile=host.pathXYZ,
                    moleculeInfo=moleculeInfo,
                    parallelize=nprocs,
                    docker=docker,
                    geom=geom)
    
    # Metadata returned to facilitate further docking, reusing products of earlier docking as hosts
    return inp_file_path, calculate_charge([guest.charge, host.charge]), calculate_multiplicity([guest.multiplicity, host.multiplicity])

def process_docking_output(inp_file_path, curr_charge, curr_multiplicity, guest, host, restraints_abs, course_key):
    course_name = course_key[0]
    # The following function splits mutli-xyz output of docker,
    # evaluates each result for satisfaction of constraints,
    # and writes only the correct outputs to separate single-xyz files
    ok_results = extract_ok_docker_results(f"{inp_file_path}.docker.struc1.all.optimized.xyz", host.n_atoms, restraints_abs, logger=logging)
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

        for col in ["ATOM", "ATOM_NAME", "RES_NAME", "CHAIN_ID", "RES_ID", "OCCUPANCY", "BETAFACTOR", "ELEMENT", "FLAVOUR", "ING"]:
            df.loc[:n_host - 1, col] = host.df[col].values
            df.loc[n_host:n_host + n_guest - 1, col] = guest.df[col].values
        # reindex residues and atoms in guest to continue the same chain
        df.loc[0:n_host + n_guest - 1, "ATOM_ID"] = range(1, n_host+n_guest+1)
        df.loc[n_host:n_host + n_guest - 1, "RES_ID"] = host.df["RES_ID"].max() + 1
        # transfer dish labels
        df.loc[:n_host - 1, "DISH"] = host.df["DISH"].values
        df.loc[n_host:n_host + n_guest - 1, "DISH"] = course_name
        # cast column types and order
        df = df[col_order].astype(col_types)
        logging.verbose( f"Result of this docking step:\n{df}")

        # Save PDB
        df_to_save = df.loc[:, ~df.columns.isin(["FLAVOUR", "ING", "DISH"])]
        new_pathPDB = new_pathXYZ.replace(".xyz", ".pdb")
        df2pdb(df=df_to_save, outPDB=new_pathPDB, remarks=[f"Eopt= {eopt} (Eh), Einter= {einter} (kcal/mol)"], logger=logging)
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
            df=df,
            restraints=restraints_abs
        )
        results_map[new_course_key] = product_obj
        c+=1
    return results_map

def optimise(ing, orca, course_name):

    # TODO separate dock and opt orca settings in config?
    orcapath, _, qmMethod_opt, _, _, _, _, _, nprocs = define_orca_params(orca)

    workdir = os.path.dirname(ing.pathXYZ)

    # Prepare ORCA input
    inp_file_path = os.path.join(workdir, "opt")
    opt_pathXYZ = f"{inp_file_path}.xyz"
    write_geom_opt_input(inp_file_path, ing, qmMethod_opt, nprocs, course_name)

    logging.info(f"Optimisation input written: {inp_file_path}.inp")
    result = run_orca(inp_file_path, orcapath, timeout=180) # TODO Put timeout in config
    logging.info(f"Optimisation complete: {inp_file_path}.out")
    if result == 25: # issue with internal coordinates breaking when angle approaches linear
        # check if opt.xyz exists - if yes, use, and maybe even proceed
        logging.info("Optimisation failed due to issues with ORCA internal coordinates as the angle approaches linearity. Evaluating outputs for usability.")
        if os.path.exists(opt_pathXYZ):
            opt_ing = process_optimisation_output(opt_pathXYZ, ing)
            if opt_ing:
                return opt_ing
        
        # Rerun with COPT (in cartesian coordinates) - but then can't restrain anything
        logging.info(f"Retrying optimisation in Cartesian coordinates")
        inp_file_path = os.path.join(workdir, "copt")
        title = f"ORCA Geometry Optimisation\n # of {ing.name}"
        moleculeInfo = {"charge": ing.charge, "multiplicity": ing.multiplicity}
        make_orca_input(orcaInput=inp_file_path,
            title=title,
            simpleInputLine=[qmMethod_opt, "COPT"],
            inputFormat="xyzfile",
            inputFile=ing.pathXYZ,
            moleculeInfo=moleculeInfo,
            parallelize=nprocs,
            geom=None)
        logging.info(f"Cartesian optimisation input written: {inp_file_path}.inp")
        result = run_orca(inp_file_path, orcapath, timeout=300) # TODO Put timeout in config
        logging.info(f"Cartesian optimisation complete: {inp_file_path}.out")
        if result == 0:
            opt_pathXYZ = f"{inp_file_path}.xyz"
            if os.path.exists(opt_pathXYZ):
                opt_ing = process_optimisation_output(opt_pathXYZ, ing)
                if opt_ing:
                    return opt_ing
        else:
            logging.error(f"Optimisation returned status code {result}. Skip to next theozyme...")
            return None

    elif result != 0 :
        logging.error(f"Optimisation returned status code {result}. Skip to next theozyme...")
        return None
    else:
        return process_optimisation_output(opt_pathXYZ, ing)

def write_geom_opt_input(inp_file_path, ing, qmMethod, nprocs, course_name):
    title = f"ORCA Geometry Optimisation\n # of {ing.name}"
    moleculeInfo = {"charge": ing.charge, "multiplicity": ing.multiplicity}

    restraints_abs = ing.restraints # TODO what happens if there's no restraints?
    curr_restraints = [r for r in restraints_abs if r.step == course_name]
    prev_restraints = [r for r in restraints_abs if r.step != course_name]

    geom_keep = []
    geom_scan = []
    geom = {}

    coords = ing.df[["X", "Y", "Z"]].to_numpy()

    # Previous restraints: simple KEEP
    for rp in prev_restraints:
        atoms = [int(s.idx) for s in rp.sele]
        if rp.property == "angle":
            val = evaluate_angle(coords, atoms[0], atoms[1], atoms[2])
        else:
            val = rp.params.val
        geom_keep.append({
            "atoms": atoms,
            "val": val
        })

    # Current restraints: SCAN
    for rc in curr_restraints:
        atoms = [int(s.idx) for s in rc.sele]

        if rc.property == "distance":
            start_val = evaluate_distance(coords, atoms[0], atoms[1])
        elif rc.property == "angle":
            start_val = evaluate_angle(coords, atoms[0], atoms[1], atoms[2])
        else:
            raise ValueError(f"Unknown restraint property: {rc.property}")

        geom_scan.append({
            "atoms": atoms,
            "start": start_val,
            "end": rc.params.val,
            "iter": max(1, round(abs(start_val - rc.params.val) / 10))
        })

    if geom_keep:
        geom["keep"] = geom_keep
    if geom_scan:
        geom["scan"] = geom_scan

    make_orca_input(orcaInput=inp_file_path,
                    title=title,
                    simpleInputLine=[qmMethod, "Opt"],
                    inputFormat="xyzfile",
                    inputFile=ing.pathXYZ,
                    moleculeInfo=moleculeInfo,
                    parallelize=nprocs,
                    geom=geom)
    
def process_optimisation_output(pathXYZ, ing):
    # Load XYZ → validate restraints → update ing → write PDB → return ing
    pathPDB = pathXYZ.replace(".xyz", ".pdb")

    if not os.path.exists(pathXYZ):
        return None

    _, opt_comment, opt_df = xyz2df(pathXYZ, logger=logging)
    coords = opt_df[["X", "Y", "Z"]].to_numpy()

    if not evaluate_restraints(coords, ing.restraints, logger=logging):
        return None

    ing.pathXYZ = pathXYZ
    ing.df[["X", "Y", "Z"]] = coords
    ing.eopt = float(opt_comment.split()[-1])
    df2pdb(
        df=ing.df,
        outPDB=pathPDB,
        remarks=[f"Eopt= {ing.eopt} (Eh), Einter= {ing.einter} (kcal/mol)"],
        logger=logging
    )
    ing.pathPDB = pathPDB
    return ing

def run_orca(input, orcapath, timeout):
    stdout_file = Path(f"{input}.out")
    stderr_file = Path(f"{input}.err")

    try:
        with stdout_file.open("w") as out, stderr_file.open("w") as err:
            # IMPORTANT: remove check=True
            result = subprocess.run(
                [orcapath, f"{input}.inp"],
                check=False,
                timeout=timeout,
                stdout=out,
                stderr=err
            )
    except subprocess.TimeoutExpired:
        logging.error(f"Process exceeded {timeout} seconds and was terminated.")
        return 1

    except Exception as e:
        # This handles only unexpected exceptions (file IO issues, etc.)
        logging.error(f"Unexpected exception during ORCA run: {e}")
        with stdout_file.open("r", errors="ignore") as f:
            lines = f.readlines()
        for line in lines[-13:]:
            logging.error(line.rstrip())
        return 1

    return result.returncode

# ----------------------------
# RESTRAINTS
# ----------------------------

def update_restraint_values(restraints, coords):
    """Update restraint target values using the optimised coordinates."""
    for r in restraints:
        atoms = [s.idx for s in r.sele]
        if r.property == "distance":
            new_val = evaluate_distance(coords, atoms[0], atoms[1])
        elif r.property == "angle":
            new_val = evaluate_angle(coords, atoms[0], atoms[1], atoms[2])
        else:
            logging.warning(f"Unknown restraint type: {r.property}")
            continue
        r.params.val = new_val

def assign_restraint_idx(guest, host, restraint, course_name):
    def _flavour_to_idx(parent, idx_ref):
        match = [(i, row["ATOM_NAME"], row["FLAVOUR"]) for i, row in parent.df.iterrows() if idx_ref in row["FLAVOUR"]]
        if not match:
            logging.error(f"Ingredient {parent.name} has no atoms matching flavour of restraint {idx_ref} for course {course_name}\n")
        else:
            logging.verbose(f"Expanding restraints for {parent.name}:")
            logging.verbose(f"\n{parent.df}")
            logging.debug(f"Original idx to be restrained: {idx_ref}")
            logging.debug(f"Guest atoms matching the restraint (idx, ATOM_NAME, FLAVOUR): {match}")
        return match

    matches = {}
    for i, atom in enumerate(restraint.sele):
        atom_parent = atom.parent # must be guest or host - no other option for iterative docker system
        atom_idx = atom.idx # at this point it is still a string
        parent_object = guest if atom_parent == "guest" else host
        match = _flavour_to_idx(parent_object, atom_idx)
        matches[i] = {"parent": atom_parent, "match": match}

    return matches

def expand_restraint(guest, host, restraint, course_name):
    # get idx matches for each atom flavour in restraint.sele
    matches = assign_restraint_idx(guest, host, restraint, course_name)

    # extract match-lists and parent names (in order)
    match_lists = [entry["match"] for entry in matches.values()]
    parents = [entry["parent"] for entry in matches.values()]

    # If any atom has zero matches → no possible restraint
    if any(len(lst) == 0 for lst in match_lists):
        return []

    # N-way Cartesian product
    expanded = []
    for combo in itertools.product(*match_lists):
        atoms = []
        for (idx, name, flav), parent in zip(combo, parents):
            atoms.append(Selection(parent, idx))

        expanded.append({
            "property": restraint.property,
            "atoms": atoms,  # 2 for distance, 3 for angle
            "val": restraint.params.val,
            "uptol": restraint.params.uptol,
            "downtol": restraint.params.downtol,
            "force": restraint.params.force,
            "orig": restraint
        })

    return expanded

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
        # Collect expansions for ALL restraints, regardless of property
        params_per_restraint = []

        for restraint in course.restraints or []:
            if restraint is None:
                continue

            expanded = expand_restraint(guest, host, restraint, course.name)

            if not expanded:
                logging.error(
                    f"Restraint {restraint.property} for course {course.name} "
                    f"could not be satisfied for guest {guest.name} / host {host.name}"
                )
                allOk = False
                params_per_restraint = []
                break

            params_per_restraint.append(expanded)

        # If no restraints → single unrestrained Course copy
        if not params_per_restraint:
            new_course = copy.deepcopy(course)
            new_course.host = host
            new_course.guests = [guest]
            new_course.restraints = []
            course_key = (str(course.name), str(guest.name), "0")
            expanded_course[course_key] = new_course
            continue

        # Cartesian product across ALL restraints
        # If one ingredient has multiple atoms of same flavour (e.g. multiple h_donor),
        # this creates multiple combinations (a "mash-up"), one for each of the atoms
        for mash_idx, mash in enumerate(itertools.product(*params_per_restraint)):
            # 'mash' is a tuple of option tuples, one per restraint
            # Build new restraint tuples representing these concrete restraint
            mashed_restraints = []
            for expanded_restraint in mash:
                orig = expanded_restraint["orig"]
                new_restraint = copy.deepcopy(orig)
                new_restraint.sele = expanded_restraint["atoms"]
                mashed_restraints.append(new_restraint)

            # Make new course copy and set host/guest to the concrete candidates and restraints
            new_course = copy.deepcopy(course)
            new_course.host = host
            new_course.guests = [guest]
            new_course.restraints = mashed_restraints

            # Unique key: include the combination index so expansion products of same host and guest don't overwrite themselves
            course_key = (str(course.name), str(guest.name), str(mash_idx))
            expanded_course[course_key] = new_course

    for course_key, course_obj in expanded_course.items():
        logging.debug(f"========= Course {course_key} restraints:")
        for i, r in enumerate(course_obj.restraints):
            logging.debug(f"----- restr {i}")
            logging.debug(r)
    return expanded_course, allOk

# ----------------------------
# MAIN LOOP
# ----------------------------
def get_molecule_value(inner_set):
    # Find all elements of length 1 (the "value" candidates)
    values = [x for x in inner_set if len(x) == 1]
    if not values:
        # No single-character value found
        return None
    elif len(values) > 1:
        # Unexpected: more than one value-like element
        raise ValueError(f"Multiple value candidates found in {inner_set}")
    return values[0]

def is_first_step(step):
    unique_molecules = {get_molecule_value(item) for item in step}
    unique_molecules.discard(None)
    if len(unique_molecules) == 1:
        return True
    return False

def row_to_signature(row, lookup):
    return "|".join(
        "x" if site is None else lookup.get_site_token(site)
        for site in row
    )

def motif_to_signature(motif, n_cols, lookup):
    row = [None] * n_cols
    for c, v in motif:
        row[c] = v
    lookup_row = [
        lookup.get(x) if x != "x" else "x"
        for x in row
    ]
    signature = "".join(lookup_row)
    return signature

def main(configPath):
    allOk = True
    logging.info(f"Cooking begins ({os.path.abspath(configPath)})")
    # Read config file
    outdir, verbosity, rmsd_threshold, orca = parse.get_default_parameters() # TODO write parameters on export from frontend
    
    config = parse.load_config(configPath)
    
    # ID lookups
    lookUp = lookup.IDLookup()
    flavour_ids = list(config.recipes[0].keys())
    idx_to_flavour = dict(enumerate(flavour_ids))
    flavour_to_idx = {v: k for k, v in idx_to_flavour.items()}
    n_cols = len(config.flavours)
    def get_flav_ids(flavour_columns):
        return [idx_to_flavour[flav] for flav in flavour_columns]
    
    
    steps, dag, roots, rows = workflow.build_execution_plan(config.recipes)
    workflow.print_execution_steps(steps, rows)
    if len(roots) > 1:
        # TODO what to do then?
        print(f"Warning: {len(roots)} roots determined.")
    first = True
    for host, guest in steps:
        # Find host file
        if first:
            host_id, host_flav_cols = workflow.determine_host(host, guest) # TODO must check that there is only one unique
            host_files = [structure.get_molec_xyz_path(cwd, config.ingredients[host_id.molecule_id].filepath)]
            first = False
        else:
            host_row = workflow.motif_to_row(host, n_cols)
            host_signature = row_to_signature(host_row, lookUp)
            host_dir = structure.prep_assembly_dir(outdir, host_signature)
            host_files = [host_dir / "best.xyz"] # TODO in the future do list(host_dir.glob("passed_*.xyz"))
            host_flav_cols = [i for i, col in enumerate(host_row) if col != 'x']
            
        # Find guest file
        guest_id, guest_flav_cols = workflow.determine_guest(host, guest)
        guest_file = structure.get_molec_xyz_path(cwd, config.ingredients[guest_id.molecule_id].filepath)
        
        # Dock dir (for output of this assembly)
        guest_row = workflow.motif_to_row(guest, n_cols)
        print(guest_row)
        guest_signature = row_to_signature(guest_row, lookUp)
        
        assembly_output_dir = structure.prep_assembly_dir(outdir, guest_signature)

        # Find restraints to apply in this step
        host_flavours = get_flav_ids(host_flav_cols)
        guest_flavours = get_flav_ids(guest_flav_cols)
        flavours = host_flavours + guest_flavours
        relevant_restraints = [
            restr for restr in config.restraints.values()
            if all(connection in flavours for connection in restr.connections)
        ]
        # New restraints are introduced as inter-host/guest bias in orca's %docker block 
        new_restraints = [
            restr for restr in relevant_restraints
            if any(connection in host_flavours for connection in restr.connections)
            and any(connection in guest_flavours for connection in restr.connections)
        ]
        # Other restraints are written in orca's %geom block
        # guest-guest restraints need to be considered in subsequent docker steps where guest has become host
        old_restraints = [
            restr for restr in relevant_restraints
            if restr in host_flavours
        ]

        flavour_map = lookup.build_flavour_map(guest_row, idx_to_flavour)
        for restr in old_restraints:
            old_atom_pairs = [
                flavour_map[fid] for fid in restr.connections
            ]
        for restr in new_restraints:
            new_atom_pairs = [
                flavour_map[fid] for fid in restr.connections
            ]
        print(new_atom_pairs)
        
        # dock
        for host_file in host_files:
            dock(host_file, guest_file, new_restraints, old_restraints, assembly_output_dir, orca)
        sys.exit(1)

        # check against restraints
        
        # optimise
        opt_ing = optimise(ing_copy, orca, course_desc.name)

        
 
    # Sort by einter value (ascending)
    leftovers_list = sorted(leftovers, key=lambda x: x.einter)
    for ing in leftovers_list:
        leftovers_dict[ing_dish_signature(ing.df)].append(ing)

    # Create output directory
    results_dir = os.path.join(outdir, "results")
    specials_dir = os.path.join(outdir, "specials")
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(specials_dir, exist_ok=True)

    logging.info(f"Saving all successfully cooked theozymes to {results_dir}...")
    # Copy files with renamed names based on order
    for i, ing in enumerate(leftovers_list, start=1):
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
    write_multi_pdb(result_paths, os.path.join(results_dir, "merged.pdb"))

    # Reduce on RMSD
    logging.info(f"Saving diverse theozymes with RMSD >= {rmsd_threshold} to {specials_dir}...")
    prod_count = 1
    dish_count = 1
    for meal in leftovers_dict.values():
        logging.debug(f"DEBUG processing dish {dish_count}")
        dish_count += 1
        diverse = []
        for ing in meal:
            logging.debug(f"DEBUG processing ingredient {ing_dish_signature(ing.df)}")
            # Compare only to already accepted structures
            isUnique = True
            min_rmsd = float('inf')
            for ref_ing in diverse:
                rmsd_val = float(rmsd.calculate_rmsd.main([ing.pathXYZ, ref_ing.pathXYZ]))
                if rmsd_val < min_rmsd:
                    min_rmsd = rmsd_val
                if rmsd_val < rmsd_threshold:
                    logging.debug(
                        f"Not saving {ing.name}: too similar to {ref_ing.name} (RMSD={rmsd_val:.3f})"
                    )
                    isUnique = False
                    break

            # If unique relative to ALL saved structures → keep it
            if isUnique:
                diverse.append(ing)
                new_name = f"result{prod_count}"
                new_pathPDB = os.path.join(specials_dir, f"{new_name}.pdb")
                shutil.copy(ing.pathPDB, new_pathPDB)
                logging.info(
                    f"Unique theozyme found: saved {ing.name} as {new_name}.pdb (Eopt={ing.eopt}, Einter={ing.einter}, min rmsd {min_rmsd})"
                )
                prod_count += 1

    # Merge into one multi-frame PDB file for easier analysis
    specials_paths = [os.path.join(specials_dir, x) for x in os.listdir(specials_dir) if x.lower().endswith(".pdb")]
    write_multi_pdb(specials_paths, os.path.join(specials_dir, "merged.pdb"))

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
    parser.add_argument('--graph_config', type=str, default=None, help="Path to the configuration YAML file generated by graph-based visual scripting website.")
    args = parser.parse_args()
    if not args.graph_config:
        # TODO del static assignemnt when finished testing
        args.graph_config = "/home/mchrnwsk/prometheozyme/graph-config.json"
    main(args.graph_config)