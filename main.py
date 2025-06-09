import yaml
import logging
import os
from datetime import datetime
import subprocess

from utils import Indices, Ingredient, Constraint, Role, update_guest_constraints, reduce_guests, print_reduced
from docking_box import calculate_docking_box
from merge_xyzs import merge_xyz

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

logger.info("Starting main script.\n")

def setup(config):
    ingredients = config.get("ingredients", [])
    roles = config.get("roles", [])

    logger.debug(f"""Ingredients are:
                {ingredients}\n""")
    logger.debug(f"""Roles are:
                {roles}\n""")

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
                {ingredient_map}\n""")

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
                {role_objects}\n""")
    # Update constraints for each role
    updated_guests_constraints = update_guest_constraints(role_objects)
    logger.debug(f"""Updated guests are:
                {updated_guests_constraints}\n""")
    # Reduce guests to unique based on constraints
    unique_guests_constraints = reduce_guests(updated_guests_constraints)
    logger.debug(f"""Unique guests are:
                {unique_guests_constraints}\n""")
    # Print missing guests after reduction for debugging
    print_reduced(updated_guests_constraints, unique_guests_constraints, logger=logger)

    return role_objects, host, ingredient_map, unique_guests_constraints

def dock(host, ingredient, outdir):
    receptor = host.path
    ligand = ingredient.path
    outdir = os.path.join(outdir, ingredient.name)
    os.makedirs(outdir, exist_ok=True)  # Ensure output directory exists
    center_x, center_y, center_z, size_x, size_y, size_z = calculate_docking_box(receptor, ligand, padding = 5.0)
    scoring = "vina"
    cnn_scoring = None
    addH = False
    exh = 8
    num_modes = 10000
    cmd = f"/opt/gnina -r {receptor} -l {ligand} --center_x {center_x} --center_y {center_y} --center_z {center_z} --size_x {size_x} --size_y {size_y} --size_z {size_z} --scoring {scoring} --cnn_scoring {cnn_scoring} --pose_sort_order Energy -o {outdir}/out.xyz --atom_terms {outdir}/atom_terms --addH {addH} --exhaustiveness {exh} --num_modes {num_modes} --quiet"
    subprocess.run(cmd.split(), check=True)
    logger.info(f"Docking for {ingredient.name} completed. Results saved in {outdir}\n")
    return f"{outdir}/out.xyz"  # Return the path to the docked output file

def main(configPath):
    configPath = "/home/mchrnwsk/theozymes/config.yaml"
    # Load the YAML file
    with open(configPath, "r") as file:
        config = yaml.safe_load(file)
    logger.info(f"Loaded configuration from {configPath}\n")
    paths = config.get("paths", {})
    workdir = paths.get("workdir", ".")
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    outdir = os.path.join(workdir, f"output_{timestamp}")

    # Prepare ingredients and roles from the configuration    
    roles, host, ingredient_map, unique_guests_constraints = setup(config)

    # Run gnina for each ingredients (constraints not evaluated here)
    for ingredient in ingredient_map.values():
        if ingredient.name == 'host':
            continue
        logger.info(f"Running gnina for ingredient: {ingredient.name}")
        docked_output = dock(host, ingredient, outdir)
        merge_xyz(host.path, docked_output, os.path.join(outdir, f"docked_{ingredient.name}.xyz"))
    logger.info(f"Docking completed. Results saved in {outdir}\n")

    # Evaluate constraints for each ingredient (list of unique guests)

    # Satisfy roles according to their priorities

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Process some ingredients and roles.")
    parser.add_argument('--config', type=str, default='config.yaml', help='Path to the configuration YAML file')
    args = parser.parse_args()
    main(args.config)