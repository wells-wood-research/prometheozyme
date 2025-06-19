import os
import re
import glob
import numpy as np
import logging
import yaml # for setup_config and setup_ingredients
import sys # To ensure custom classes are found if not run as a script

from utils import read_xyz, get_atom_count, convert_optimised_arr_xyz_to_pdb
from main import setup_config, setup_ingredients
from glycinate import add_glycine_to_pdb

# Configure logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

# Console handler for logging output
ch = logging.StreamHandler(sys.stdout)
ch.setLevel(logging.INFO) # Set to INFO for general messages, DEBUG for detailed
ch.setFormatter(formatter)
logger.addHandler(ch)

def process_arrangement_directories(base_directory, config_path):
    """
    Finds arrangement files, extracts guest info, and converts pull.xyz to PDB.

    Args:
        base_directory (str): The base directory to search in.
        config_path (str): Path to the configuration YAML file for ingredient definitions.
    """
    logger.info(f"Starting processing in directory: {base_directory}")

    # Call setup_config from imported main module
    config, _, _, _, _, _ = setup_config(config_path)
    if not config:
        logger.error("Failed to load configuration. Aborting PDB conversion.")
        return

    # Call setup_ingredients from imported main module
    _, host, ingredient_map, _ = setup_ingredients(config)
    if not host or not ingredient_map:
        logger.error("Failed to setup ingredients from config. Aborting PDB conversion.")
        return
    
    # Call get_atom_count from imported utils module
    host_atom_count = get_atom_count(host.path, logger=logger)
    if host_atom_count == 0:
        logger.error(f"Could not get atom count for host from '{host.path}'. Please check host file validity. Aborting.")
        return

    # Find all arrangement_*.xyz files using glob for pattern matching
    # The pattern 'arrangement_[0-9]*.xyz' ensures we match 'arrangement_1.xyz', 'arrangement_123.xyz' etc.
    arrangement_files = glob.glob(os.path.join(base_directory, "arrangement_[0-9]*.xyz"))
    
    if not arrangement_files:
        logger.info(f"No 'arrangement_*.xyz' files found in '{base_directory}'.")
        return

    logger.debug(f"Found arrangement_*.xyz files: {[os.path.basename(f) for f in arrangement_files]}")

    for arr_file in sorted(arrangement_files): # Process in sorted order for predictability
        # Extract arrangement number (e.g., '1' from 'arrangement_1.xyz') and regex to ensure it's a digit
        match = re.search(r'arrangement_(\d+)\.xyz$', os.path.basename(arr_file))
        if not match:
            logger.warning(f"Skipping '{os.path.basename(arr_file)}': does not match expected 'arrangement_N.xyz' pattern.")
            continue
        arrangement_num = match.group(1)

        # Construct the path to the corresponding pull/pull.xyz file
        pull_xyz_path = os.path.join(base_directory, f"arrangement_{arrangement_num}", "pull", "pull.xyz")

        if not os.path.exists(pull_xyz_path):
            logger.warning(f"Skipping arrangement {arrangement_num}: Corresponding pull file '{pull_xyz_path}' not found.")
            continue

        logger.info(f"\n--- Processing Arrangement {arrangement_num} ---")
        logger.info(f"  Reading Arrangement XYZ: {os.path.basename(arr_file)}")
        logger.info(f"  Associated Pull XYZ: {os.path.relpath(pull_xyz_path, base_directory)}")

        # Read the comment line from arrangement_*.xyz to get guest order
        # (expect only one structure in arrangement_N.xyz)
        _, arr_comment, _, _ = read_xyz(arr_file, logger=logger)[0]
        if not arr_comment:
            logger.warning(f"Could not read data from '{arr_file}'. Skipping this arrangement.")
            continue

        # Extract guest names from the comment line.
        # The regex captures the guest name (e.g., "ser", "arg") which is the first group (\w+).
        # re.findall with this pattern will return a list of tuples like [('ser', 'base@0'), ('arg', 'stabSubO@45')]
        guest_matches = re.findall(r'(\w+)\(([^)]+)\@\d+\)', arr_comment)
        
        # Extract just the guest name (first element of each tuple)
        guest_names_ordered = [match[0] for match in guest_matches]
        
        if not guest_names_ordered:
            logger.info(f"No specific guest names found in comment for {os.path.basename(arr_file)}: '{arr_comment}'. "
                        f"PDB will only contain host and un-named remaining atoms if any.")
            # If no guests found in comment the final PDB would be only the host
            continue

        logger.info(f"  Extracted guest names (in order): {guest_names_ordered}")

        # Read the pull.xyz file, which contains the combined host and guests
        # (again expected to contain a single combined structure)
        pull_xyz_structure = read_xyz(pull_xyz_path, logger=logger)[0]
        if not pull_xyz_structure:
            logger.warning(f"Could not read data from '{pull_xyz_path}'. Skipping PDB conversion for this arrangement.")
            continue

        # Define output PDB path
        output_pdb_dir = os.path.join(base_directory, "final_theozymes", "no_gly")
        output_gly_dir = os.path.join(base_directory, "final_theozymes", "gly")
        os.makedirs(output_pdb_dir, exist_ok=True) # Ensure the output directory exists
        os.makedirs(output_gly_dir, exist_ok=True) # Ensure the output directory exists
        output_pdb_path = os.path.join(output_pdb_dir, f"arrangement_{arrangement_num}.pdb")
        output_gly_path = os.path.join(output_gly_dir, f"arrangement_{arrangement_num}.pdb")

        # Call the function to write the PDB file
        convert_optimised_arr_xyz_to_pdb(output_pdb_path, pull_xyz_structure, host_atom_count, guest_names_ordered, ingredient_map, logger=logger)
        add_glycine_to_pdb(output_pdb_path, output_gly_path, glycine_pdb="/home/mchrnwsk/theozymes/ingredients/glycine.pdb")

    logger.info("\nAll arrangements processed.")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Convert arrangement_N/pull/pull.xyz files to PDB with correctly assigned residue names."
    )
    parser.add_argument(
        '--dir', 
        type=str,
        help="Path to the base output directory (e.g., /home/mchrnwsk/theozymes/docking/output_2025-06-17_12-00-00/)"
    )
    parser.add_argument(
        '--config', 
        type=str, 
        default='config.yaml', # Default value assumes config.yaml is in the same directory as this script
        help="Path to the configuration YAML file (e.g., /path/to/config.yaml). "
             "Defaults to 'config.yaml' in the script's directory."
    )

    args = parser.parse_args()
    dir = args.dir
    if not dir:
        dir = "/home/mchrnwsk/theozymes/docking/output_2025-06-17_12-00-00"
    if not os.path.exists(dir):
        logger.error(f"Dir {dir} does not exist.")

    config = args.config
    if not config:
        config = "/home/mchrnwsk/theozymes/config.yaml"
    if not os.path.exists(config):
        logger.error(f"Config {config} does not exist.")

    # Execute the main processing function
    process_arrangement_directories(dir, config)
