import yaml
import logging
import os
from datetime import datetime
import subprocess
import shutil

from define import Indices, Ingredient, Constraint, Role, update_guest_constraints, reduce_guests, print_reduced
from dock import dock
from utils import read_xyz, merge_xyz, split_multi_pdb, write_multi_pdb, append_scores, get_atom_count
from evaluate import filter_conformations
from arrange import arrange_guests
from protonate import protonate_all, deprotonate_selected
from reindex import reindex
from optimise import optimise

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
    configPath = "/home/mchrnwsk/theozymes/config.yaml"
    # Load the YAML file
    with open(configPath, "r") as file:
        config = yaml.safe_load(file)
    logger.info(f"Loaded configuration from {configPath}\n")
    misc = config.get("misc", {})
    dock_params = config.get("docking", {})
    redock_params = config.get("redocking", {})
    backbone = misc.get("backbone", False)
    evaluate_backbone = misc.get("evaluate_backbone", False)
    workdir = misc.get("workdir", ".")
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    outdir = os.path.join(workdir, f"output_{timestamp}")

    return config, outdir, dock_params, redock_params, backbone, evaluate_backbone

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

def process_reindexing(protonated_structures, reindex_reference_path, reindex_output_dir, ing_name, ing_id, logger):
    successful = []
    failed = []
    for idx, path in enumerate(protonated_structures):
        reindexed_pdb_path = reindex(reindex_reference_path, path, reindex_output_dir, logger=logger)
        if reindexed_pdb_path is None:
            if logger:
                logger.info(f"Reindexing failed for pose {idx} of ingredient {ing_name}")
            continue
        successful.append(reindexed_pdb_path)

    # Check if all poses failed
    if not successful:
        failed.append(ing_id)
        if logger:
            logger.warning(f"No poses were successfully reindexed for ingredient {ing_name}")
    else:
        if logger:
            logger.info(f"For ingredient {ing_name}, {len(successful)} out of {len(protonated_structures)} poses were successfully reindexed")

    return successful, failed

def process_docking(ing, host, outdir, dock_params, redocking=False):
    # Handles the docking, per-conformation protonation, reindexing, and deprotonation for a single ingredient.
    logger.info(f"Processing docking for {ing.name}{' (redocking)' if redocking else ''}...")
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
        logger.warning(f"Docked output for {ing.name} not found or is empty at {docked_path}. Skipping further processing.")
        return None
    
    # Protonate output of initial docking on all atoms
    prot_path, prot_atom_count = protonate_all(docked_path)
    if not os.path.exists(prot_path) or os.path.getsize(prot_path) == 0:
        logger.warning(f"Protonated output for {ing.name} not found or is empty at {prot_path}. Skipping further processing.")
        return None

    # Prepare reference for reindexing (PDB version of original ingredient)
    reindex_reference_path = f"{os.path.splitext(ing.path)[0]}_ref.pdb"
    if not os.path.exists(reindex_reference_path):
        logger.warning("Reference ingredient in PDB format not found - using XYZ instead; consider adding a PDB reference with same atom indices and path as XYZ.")
        reindex_reference_path = ing.path

    # Split the protonated multi-PDB into individual PDB files
    protonated_structures = split_multi_pdb(prot_path, prot_output_dir, logger)
    # Reindex protonated PDB files individually
    successful, failed = process_reindexing(protonated_structures, reindex_reference_path, reindex_output_dir, ing.name, ing.id, logger)
    if failed:
        logger.error(f"Reindexing failed for {ing.name}. Investigate output manually.")
        logger.warning(f"Failed ingredients: {', '.join(failed)}")
        return None, failed
    # Merge protonated and reindexed PDB files into one multi-PDB file
    multi_prot_path = os.path.join(outdir, ing.name, f"{ing.name}_docked_prot.pdb")
    write_multi_pdb(successful, multi_prot_path)

    # Deprotonate specific atoms
    deprot_path, deprot_atom_count = deprotonate_selected(multi_prot_path, getattr(ing.indices, "deprotonate", []))
    if not os.path.exists(deprot_path):
        logger.error(f"Deprotonated file not found at expected path: {deprot_path}.")
    if deprot_atom_count != reference_atom_count:
        failed.append(ing.name)
        logger.error(f"Deprotonation failed for {ing.name}. Investigate output manually.")
        logger.warning(f"Failed ingredients: {', '.join(failed)}")
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
    logger.info(f"Final merged host-guest file saved to: {final_merged_host_guest_path}")

    # Append scores from gnina docking
    scores_file_from_dock = os.path.join(os.path.dirname(docked_path), "scores.txt")
    if os.path.exists(scores_file_from_dock):
        append_scores(final_merged_host_guest_path, scores_file_from_dock, logger)
        logger.info(f"Appended scores to {final_merged_host_guest_path}")
    else:
        logger.warning(f"Scores file not found at {scores_file_from_dock}. Skipping score appending.")

    # Clean up temporary processing directory
    shutil.rmtree(temp_processing_dir)
    logger.debug(f"Cleaned up temporary directory: {temp_processing_dir}")

    return final_merged_host_guest_path, failed

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
        redocked_path, failed = process_docking(ing, host, outdir, redock_dock_params, redocking=True)

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
    config, outdir, dock_params, redock_params, backbone, evaluate_backbone = setup_config(configPath)

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
    failed = []
    for ing in ingredient_map.values():
        if ing.name == 'host':
            continue
        logger.info(f"Running gnina for ingredient: {ing.name}")
        outpath, failed_i = process_docking(ing, host, outdir, dock_params)
        if failed_i:
            failed.append(failed_i)
    logger.info(f"Docking completed. Results saved in {outdir}\n")
    
    # Evaluate constraints for each ingredient (list of unique guests)
    # for ingredient in unique_guests_constraints, delete in copy conformations that don't satisfy constraints
    for ing in unique_guests_constraints:
        logger.info(f"\nEvaluating constraints for ingredient: {ing.name}, role: {ing.role_title}")

        docked_path = os.path.join(outdir, f"{ing.name}.xyz")
        if not os.path.exists(docked_path):
            logger.error(f"Docked output for {ing.name} not found at {docked_path}. Skipping constraint evaluation.")
            continue

        # Filter conformations that don't satisfy constraints
        process_constraints(docked_path, ing, host, outdir, evaluate_backbone, dock_params, redock_params, logger)
        
    # Satisfy roles according to their priorities
    sorted_final_arrangements = arrange_guests(roles, unique_guests_constraints, host_atom_count, host_coords, host_atom_types, outdir, logger)

    for arr in sorted_final_arrangements:
        failed_ids = [guest["obj"] for guest in arr["guests_info"] if guest["obj"] in failed]

        if failed_ids:
            logger.warning(f"Arrangement at path {arr['path']} includes failed guests with IDs: {failed_ids}. Protonate/reindex/deprotonate manually before proceeding.")
        else:
            logger.debug(f"Arrangement at path {arr['path']} has no failed guests.")
        # simple constrained optimisation
        arr_optimised, arr_reactant = optimise(arr, host_atom_count, ingredient_map, backbone, logger)
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