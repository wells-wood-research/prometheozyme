import yaml
import logging
import os
from datetime import datetime
import subprocess

from define import Indices, Ingredient, Constraint, Role, update_guest_constraints, reduce_guests, print_reduced
from dock import dock
from utils import read_xyz, merge_xyz, append_scores
from evaluate import filter_conformations
from arrange import arrange_guests
from protonate import protonate
from reindex import reindex

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

def setup_config(configPath):
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

    return config, outdir, dock_params, redock_params, evaluate_backbone

def setup_ingredients(config):
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

def process_docking(ing, host, outdir, dock_params, redocking=False, logger=logger):
    # dock
    docked_path = dock(host, ing, outdir, dock_params, redocking, logger)

    # protonate
    protonated_path = protonate(docked_path, getattr(ing.indices, "deprotonate", [0]))

    # reindex
    reindexed_path = reindex(ing.path, protonated_path, os.path.join(outdir, ing.name))

    # merge docked protonated output with host
    merged_path = os.path.join(outdir, f"docked_{ing.name}.xyz")
    merge_xyz(host.path, protonated_path, merged_path)

    # write scores as comment line for each structure in multi XYZ merged output
    append_scores(merged_path, docked_path.replace("out.xyz", "scores.txt"), logger)

    return merged_path

def process_constraints(docked_path, ing, host, outdir, evaluate_backbone, dock_params, redock_params, logger):
    # Filter conformations that don't satisfy constraints - only valid written to filtered_path
    valid_structures, filtered_path = filter_conformations(docked_path, host.path, ing.name, ing.role_title, ing.indices, ing.constraints, evaluate_backbone, logger)

    # Repeat docking if no valid conformations found
    if len(valid_structures) == 0:
        logger.warning(f"No poses that satisfy constraints found on first docking attempt. Redocking with {redock_params}...")
        # Update dock_params with redocking parameters
        redock_dock_params = dock_params.copy()
        redock_dock_params.update(redock_params)
        # Repeat docking with redocking parameters
        redocked_path = process_docking(ing, host, outdir, redock_params, redocking=True)

        # Repeat filtering conformations that don't satisfy constraints
        valid_structures, filtered_path = filter_conformations(redocked_path, host.path, ing.name, ing.role_title, ing.indices, ing.constraints, evaluate_backbone, logger)
        if len(valid_structures) == 0:
            logger.error(f"No poses that satisfy constraints {ing.constraints} for {ing.name}, role: {ing.role_title} found on repeat docking.")
            return
    ing.update_conformations(filtered_path)
    logger.info(f"Filtered conformations are saved in {filtered_path}")
    return filtered_path
    
def main(configPath):
    # Read config file
    config, outdir, dock_params, redock_params, evaluate_backbone = setup_config(configPath)

    # Prepare ingredients and roles from the configuration    
    roles, host, ingredient_map, unique_guests_constraints = setup_ingredients(config)
    # Sort roles by priority (lower number = higher priority, e.g., -1 is highest)
    roles = sorted(roles, key=lambda x: x.priority)

    # Read key host data for reuse across application
    host_data = read_xyz(host.path, logger)
    if not host_data:
        logger.error(f"Could not read host file: {host.path}. Exiting arrangement process.")
        return []
    host_atom_count, _, host_coords, host_atom_types = host_data[0]

    # Run gnina for each ingredients (constraints not evaluated here)
    for ing in ingredient_map.values():
        if ing.name == 'host':
            continue
        logger.info(f"Running gnina for ingredient: {ing.name}")
        outpath = process_docking(ing, host, outdir, dock_params)
    logger.info(f"Docking completed. Results saved in {outdir}\n")
    
    # Evaluate constraints for each ingredient (list of unique guests)
    # for ingredient in unique_guests_constraints, delete in copy conformations that don't satisfy constraints
    for ing in unique_guests_constraints:
        logger.info(f"\nEvaluating constraints for ingredient: {ing.name}, role: {ing.role_title}")

        docked_path = os.path.join(outdir, f"docked_{ing.name}.xyz")
        if not os.path.exists(docked_path):
            logger.error(f"Docked output for {ing.name} not found at {docked_path}. Skipping constraint evaluation.")
            continue

        # Filter conformations that don't satisfy constraints
        process_constraints(docked_path, ing, host, outdir, evaluate_backbone, dock_params, redock_params, logger)
        
    # Satisfy roles according to their priorities
    sorted_final_arrangements = arrange_guests(roles, unique_guests_constraints, host_atom_count, host_coords, host_atom_types, outdir, logger)

    #for arr in sorted_final_arrangements:
        # simple constrained optimisation
        #arr_optimised = optimise(arr, arr)
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