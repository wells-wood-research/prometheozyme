import yaml
import logging
import os
from datetime import datetime
import subprocess

from define import Indices, Ingredient, Constraint, Role, update_guest_constraints, reduce_guests, print_reduced
from docking_box import calculate_docking_box
from utils import read_xyz, merge_xyz, append_scores
from evaluate import filter_conformations
from arrange import arrange_guests
from drOrca import make_orca_input

# Configure logging to output to both console and file
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

# Console handler
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(formatter)
logger.addHandler(ch)

# File handler
fh = logging.FileHandler('main.log', mode='w')
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)
logger.addHandler(fh)

def setup(config):
    ingredients = config.get("ingredients", [])
    roles = config.get("roles", [])

    logger.debug(f"""Ingredients are:
                {ingredients}""")
    logger.debug(f"""Roles are:
                {roles}""")

    ingredient_objects = []
    ingredient_map = {}  # Map ingredient names to objects for role processing
    for ing in ingredients:
        indices_dict = ing.get('indices', {})
        indices_obj = Indices(**indices_dict)  # Create Indices object with all YAML indices
        ingredient_obj = Ingredient(
            path=ing['path'],
            charge=ing['charge'],
            multiplicity=ing['multiplicity'],
            indices=indices_obj,
            name=ing['name']
        )
        ingredient_objects.append(ingredient_obj)
        ingredient_map[ing['name']] = ingredient_obj
    logger.debug(f"""Ingredient map is:
                {ingredient_map}""")

    # Find the substrate (host) ingredient, assumed to be 'sub'
    host = next((ing for ing in ingredient_objects if ing.name == 'host'), None)
    if not host:
        print("Error: No 'sub' ingredient found for role host")
        exit(1)

    # Create Role objects from YAML
    role_objects = []
    for role in roles:
        # Map candidates to Ingredient objects
        candidates = [ingredient_map[cand['name']] for cand in role['candidates']]
        
        # Handle constraints if present
        constraints_data = role.get('constraints', [])
        constraints = []
        for cons in constraints_data:
            constraint = Constraint(
                guestIdx=cons['guestIdx'],
                guestType=cons.get('guestType', 'iter'),
                hostIdx=cons['hostIdx'],
                hostType=cons.get('hostType', 'iter'),
                val=cons['val']
            )
            constraints.append(constraint)
        
        # Create Role object
        role_obj = Role(
            title=role['name'],
            priority=role['priority'],
            guests=candidates,
            host=host,
            constraints=constraints
        )
        role_objects.append(role_obj) 
    logger.debug(f"""Role objects are:
                {role_objects}""")
    # Update constraints for each role
    updated_guests_constraints = update_guest_constraints(role_objects)
    logger.debug(f"""Updated guests are:
                {updated_guests_constraints}""")
    # Reduce guests to unique based on constraints
    unique_guests_constraints = reduce_guests(updated_guests_constraints)
    logger.debug(f"""Unique guests are:
                {unique_guests_constraints}\n""")
    # Print missing guests after reduction for debugging
    print_reduced(updated_guests_constraints, unique_guests_constraints, logger)

    return role_objects, host, ingredient_map, unique_guests_constraints

def dock(host, ingredient, outdir, dock_params, redocking=False):
    receptor = host.path
    ligand = ingredient.path
    outdir = os.path.join(outdir, ingredient.name)
    os.makedirs(outdir, exist_ok=True)  # Ensure output directory exists
    
    # Calculate docking box
    center_x, center_y, center_z, size_x, size_y, size_z = calculate_docking_box(receptor, ligand, padding=5.0)
    
    # Update dock_params with docking box coordinates
    dock_params = dock_params.copy()  # Avoid modifying the original
    dock_params.update({
        'center_x': center_x,
        'center_y': center_y,
        'center_z': center_z,
        'size_x': size_x,
        'size_y': size_y,
        'size_z': size_z
    })

    # Construct the command as a list
    cmd = [
        "/opt/gnina",
        "-r", receptor,
        "-l", ligand,
        "--center_x", str(dock_params['center_x']),
        "--center_y", str(dock_params['center_y']),
        "--center_z", str(dock_params['center_z']),
        "--size_x", str(dock_params['size_x']),
        "--size_y", str(dock_params['size_y']),
        "--size_z", str(dock_params['size_z']),
        "--scoring", dock_params['scoring'],
        "--cnn_scoring", dock_params['cnn_scoring'],
        "--pose_sort_order", dock_params['pose_sort_order'],
        "-o", os.path.join(outdir, "out.xyz"),
        "--atom_terms", os.path.join(outdir, "atom_terms"),
        "--exhaustiveness", str(dock_params['exhaustiveness']),
        "--num_modes", str(dock_params['num_modes'])
    ]

    # Include optional parameters if not None or False
    if dock_params['addH']:
        cmd.append("--addH")
    if dock_params['stripH']:
        cmd.append("--stripH")
    if dock_params.get("no_gpu", False):
       cmd.append("--no_gpu")
    if dock_params.get("min_rmsd_filter", None):
        cmd.extend(["--min_rmsd_filter", str(dock_params['min_rmsd_filter'])])
    if dock_params.get("temperature", None):
        cmd.extend(["--temperature", str(dock_params['temperature'])])
    if dock_params.get("num_mc_steps", None):
        cmd.extend(["--num_mc_steps", str(dock_params['num_mc_steps'])])

    with open(os.path.join(outdir, 'scores.txt'), 'w') as outfile:
        subprocess.run(cmd, check=True, stdout=outfile, stderr=subprocess.STDOUT)
    logger.info(f"Docking for {ingredient.name} {'(redocking)' if redocking else ''} completed. Results saved in {outdir}")
    return os.path.join(outdir, "out.xyz")  # Return the path to the docked output file


def calculate_charge(charges):
    return int(sum(charges))
    
def calculate_multiplicity(multiplicities):
    # Assume molecules are weakly or non-interacting
    return int(sum(m - 1 for m in multiplicities) / 2)

def optimise(arr, path, host_atom_count):

    arrName = os.path.splitext(os.path.basename(arr["path"]))[0]
    orcaInputDir = os.path.join(os.path.dirname(arr["path"]), arrName) 
    os.makedirs(orcaInputDir, exist_ok=True)
    orcaInput = os.path.join(orcaInputDir, "opt.inp")

    title = f"""Initial optimisation of {arrName}:\n
                {arr["desc"]}"""
    
    qmMethod = "XTB2"
    inputFormat = "xyzfile"

    moleculeInfo = {"charge": calculate_charge([guest["obj"].charge for guest in arr["guests_info"]]),
                    "multiplicity": calculate_multiplicity([guest["obj"].multiplicity for guest in arr["guests_info"]])}

    parallelize = 8
    maxcore = 2500

    geom = {"keep": []}
    for guest_info in arr["guests_info"]:
        guest = guest_info["obj"]
        constraints = guest.constraints
        for constraint in guest.constraints:
            guestIdx, hostIdx, val = constraint
            guestIdx = host_atom_count + guestIdx
            atoms = [guestIdx, hostIdx]
            geom["keep"].append({"atoms": atoms, "val": val})

    make_orca_input(orcaInput = orcaInput,
                    title = title,
                    qmMethod = qmMethod,
                    inputFormat = inputFormat,
                    inputFile = path,
                    moleculeInfo = moleculeInfo,
                    parallelize = parallelize,
                    maxcore = maxcore,
                    qmmm = None,
                    geom = geom,
                    neb = None,
                    scf = None,
                    docker = None)
    arr_optimised = orcaInput.replace(".inp", ".xyz")
    return arr_optimised

def protonate(arr, path):
    return None
 
#def pull_backbone_out(arrangement, arr_optimised):
#    return arr_reactant
#
#def prepare_product(arrangement, arr_reactant):
#    return arr_product

def main(configPath):
    configPath = "/home/mchrnwsk/theozymes/config.yaml"
    # Load the YAML file
    with open(configPath, "r") as file:
        config = yaml.safe_load(file)
    logger.info(f"Loaded configuration from {configPath}\n")
    misc = config.get("misc", {})
    dock_params = config.get("docking", {})
    redock_params = config.get("redocking", {})
    evaluate_backbone = misc.get("evaluate_backbone", False)
    workdir = misc.get("workdir", ".")
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    outdir = os.path.join(workdir, f"output_{timestamp}")

    # Prepare ingredients and roles from the configuration    
    roles, host, ingredient_map, unique_guests_constraints = setup(config)
    # Sort roles by priority (lower number = higher priority, e.g., -1 is highest)
    roles = sorted(roles, key=lambda x: x.priority)

    host_data = read_xyz(host.path, logger)
    if not host_data:
        logger.error(f"Could not read host file: {host.path}. Exiting arrangement process.")
        return []
    host_atom_count, host_comment, host_coords, host_atom_types = host_data[0]

    # Run gnina for each ingredients (constraints not evaluated here)
    docked_guests = []
    for ingredient in ingredient_map.values():
        if ingredient.name == 'host':
            continue
        logger.info(f"Running gnina for ingredient: {ingredient.name}")
        docked_output = dock(host, ingredient, outdir, dock_params, redocking=False)
        merged_path = os.path.join(outdir, f"docked_{ingredient.name}.xyz")
        merge_xyz(host.path, docked_output, merged_path)
        append_scores(merged_path, docked_output.replace("out.xyz", "scores.txt"), logger)
        docked_guests.append(merged_path)
    logger.info(f"Docking completed. Results saved in {outdir}\n")
    
    # Evaluate constraints for each ingredient (list of unique guests)
    # for ingredient in unique_guests_constraints, delete in copy conformations that don't satisfy constraints
    for ing in unique_guests_constraints:
        logger.info(f"\nEvaluating constraints for ingredient: {ing.name}, role: {ing.role_title}")

        merged_path = os.path.join(outdir, f"docked_{ing.name}.xyz")
        if not os.path.exists(merged_path):
            logger.error(f"Docked output for {ing.name} not found at {merged_path}. Skipping constraint evaluation.")
            continue

        # Filter conformations that don't satisfy constraints - only valid written to filtered_path
        valid_structures, filtered_path = filter_conformations(merged_path, host.path, ing.name, ing.role_title, ing.indices, ing.constraints, evaluate_backbone, logger)
        if len(valid_structures) == 0:
            # Repeat docking with redocking parameters
            logger.warning(f"No poses that satisfy constraints found on first docking attempt. More granular redocking...")
            # Update dock_params with redocking parameters
            redock_dock_params = dock_params.copy()
            redock_dock_params.update(redock_params)
            docked_output = dock(host, ingredient_map[ing.name], outdir, redock_dock_params, redocking=True)
            merged_path = os.path.join(outdir, f"docked_{ing.name}.xyz")
            merge_xyz(host.path, docked_output, merged_path)
            valid_structures, filtered_path = filter_conformations(merged_path, host.path, ing.name, ing.role_title, ing.indices, ing.constraints, evaluate_backbone, logger)
            if len(valid_structures) == 0:
                logger.error(f"No poses that satisfy constraints {ing.constraints} for {ing.name}, role: {ing.role_title} found on repeat docking.")
                continue
        ing.update_conformations(filtered_path)
        logger.info(f"Filtered conformations are saved in {filtered_path}")

    # Satisfy roles according to their priorities
    sorted_final_arrangements = arrange_guests(roles, unique_guests_constraints, host.path, outdir, logger)

    for arr in sorted_final_arrangements:
        # protonate!
        arr_protonated = protonate(arr, arr["path"], host_atom_count)
        # simple constrained optimisation
        arr_optimised = optimise(arr, arr_protonated)
    #    # constrained optimisation with pulling backbone out
    #    arr_reactant = pull_backbone_out(arr, arr_optimised)
    #    # prepare product
    #    arr_product = prepare_product(arr, arr_reactant)

    #for arrangement in final_arrangements_list:
        # run NEB

    logger.info("""
          Main script complete!
          ⚝⭒٭⋆⚝⭒٭⋆⚝⭒٭⋆⚝⭒٭⋆⚝⭒٭⋆⚝⭒٭⋆
          """)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Process some ingredients and roles.")
    parser.add_argument('--config', type=str, default='config.yaml', help='Path to the configuration YAML file')
    args = parser.parse_args()
    main(args.config)