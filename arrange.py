import os
import numpy as np

# Assuming these are available from utils.py
from utils import get_vdw_radius, read_xyz, write_xyz, calculate_distance

# Constants for filtering and limiting output
MAX_FINAL_ARRANGEMENTS = 100  # Maximum number of unique ingredient-role combinations to keep

# Custom exception for early termination
class MaxArrangementsReached(Exception):
    """Custom exception to signal that the maximum number of arrangements has been found."""
    pass

def arrange_guests(roles, unique_guests_constraints, host_path, outdir, logger=None):
    """
    Finds a limited number of optimal and diverse arrangements of guests around a host,
    satisfying as many roles as possible according to their priority, and ensuring no
    spatial overlap between guests. For each unique combination of ingredients and roles,
    it finds and immediately outputs the arrangement with the most negative sum of conformer energies.
    The search stops when MAX_FINAL_ARRANGEMENTS unique ingredient-role combinations are found.

    Args:
        roles (list): A list of objects, each with 'name' (str), 'candidates' (list),
                      'constraints' (list), and 'priority' (int).
        unique_guests_constraints (list): A list of objects representing unique guest candidates
                                          that have valid conformations. Each object should have:
                                          'name' (str), 'role_title' (str, matching role.title),
                                          'conformations' (str, path to multi-XYZ file of valid conformations),
                                          'charge', 'multiplicity', 'id', 'indices', 'path'.
        host_path (str): Path to the host XYZ file.
        outdir (str): Output directory for saving the final arrangements.
        logger: A logger object for outputting information and warnings.

    Returns:
        list: A list of dictionaries, where each dictionary represents a valid arrangement.
              Each dictionary contains:
              - 'guests': A list of dictionaries for each selected guest, including 'name',
                          'role_title', and 'conformation_idx' (index in its multi-XYZ file).
              - 'satisfied_roles_count': The total number of roles satisfied by this arrangement.
    """

    # Group unique_guests_constraints by role_title for easier lookup
    guests_by_role_title = {}
    for guest_uc in unique_guests_constraints:
        if guest_uc.role_title not in guests_by_role_title:
            guests_by_role_title[guest_uc.role_title] = []
        guests_by_role_title[guest_uc.role_title].append(guest_uc)

    # This dictionary will store the best arrangement found so far for each unique signature.
    # Key: tuple of sorted (guest_id, role_title) pairs (the signature)
    # Value: {'total_energy_sum': float, 'guests_info': list_of_dicts_for_this_arrangement}
    best_arrangement_info_per_signature = {}
    
    max_satisfied_roles = -1
    saved_arrangements_count = 0 # Counter for naming output files sequentially

    # Read host coordinates and atom types once for combining later
    host_data = read_xyz(host_path, logger)
    if not host_data:
        logger.error(f"Could not read host file: {host_path}. Exiting arrangement process.")
        return []
    
    host_atom_count, host_comment, host_coords, host_atom_types = host_data[0]

    # Recursive helper function to explore arrangements
    def _recursively_arrange(
        current_role_index,
        current_arrangement_guests_info, # list of dicts: {'obj': guest_obj, 'coords': np.array, 'atom_types': list, 'conf_idx': int, 'energy': float}
        satisfied_roles_count,
        current_energy_sum, # Sum of energies of conformers in current_arrangement_guests_info
    ):
        nonlocal max_satisfied_roles, saved_arrangements_count

        # Base Case: All roles have been considered
        if current_role_index == len(roles):
            # If we found an arrangement that satisfies more roles, clear previous ones
            if satisfied_roles_count > max_satisfied_roles:
                max_satisfied_roles = satisfied_roles_count
                best_arrangement_info_per_signature.clear() # Clear all previous if a better count is found
                logger.info(f"Found arrangement(s) satisfying more roles ({max_satisfied_roles}). Resetting unique arrangement buffer.")

            # Only consider arrangements that satisfy the current maximum number of roles
            if satisfied_roles_count == max_satisfied_roles:
                # 1. Create the signature for the current arrangement
                # The signature identifies the unique set of guests and the roles they satisfied.
                current_signature_parts = []
                for g_info in current_arrangement_guests_info:
                    current_signature_parts.append((g_info['obj'].id, g_info['obj'].role_title))
               
                # Sort the parts to ensure a canonical representation regardless of order of addition
                current_signature = tuple(sorted(current_signature_parts)) 

                # 2. Check if this combination of ingredients/roles is new or if this is a better energy
                existing_best_info = best_arrangement_info_per_signature.get(current_signature)

                if existing_best_info is None or current_energy_sum < existing_best_info['total_energy_sum']:
                    # This is a new unique ingredient-role combination or an improved (lower energy) one
                    logger.debug(f"New/Improved arrangement found for signature {current_signature} with energy {current_energy_sum:.3f}.")
                    
                    # Store the details for output
                    arrangement_details_for_output = {
                        'guests_info': current_arrangement_guests_info,
                        'satisfied_roles_count': satisfied_roles_count,
                        'total_energy_sum': current_energy_sum
                    }
                    
                    best_arrangement_info_per_signature[current_signature] = arrangement_details_for_output

                    # Immediately output this arrangement
                    saved_arrangements_count += 1
                    arrangement_output_path = os.path.join(outdir, f"arrangement_{saved_arrangements_count}.xyz")
                    
                    # Build comment and combine coordinates for output
                    combined_coords = host_coords.tolist()
                    combined_atom_types = host_atom_types[:]
                    
                    arrangement_comment_parts = [
                        f"Arrangement {saved_arrangements_count} (Roles: {satisfied_roles_count}, Energy Sum: {current_energy_sum:.3f})"
                    ]
                    for guest_info in current_arrangement_guests_info:
                        combined_coords.extend(guest_info['coords'].tolist())
                        combined_atom_types.extend(guest_info['atom_types'])
                        arrangement_comment_parts.append(f"{guest_info['obj'].name}({guest_info['obj'].role_title}@{guest_info['conf_idx']})")

                    write_xyz(arrangement_output_path, " ".join(arrangement_comment_parts), np.array(combined_coords), combined_atom_types)
                    logger.info(f"Arrangement {saved_arrangements_count} saved to {arrangement_output_path}")

                    # Check for MAX_FINAL_ARRANGEMENTS (now refers to max unique signatures)
                    if len(best_arrangement_info_per_signature) == MAX_FINAL_ARRANGEMENTS:
                        logger.info(f"Reached MAX_FINAL_ARRANGEMENTS ({MAX_FINAL_ARRANGEMENTS}) unique signature arrangements. Stopping search early.")
                        raise MaxArrangementsReached
            return # Simply return from base case, do not control parent loop flow

        current_role = roles[current_role_index]
        
        # Option 1: Try to satisfy the current role with a candidate guest
        candidates_for_this_role = guests_by_role_title.get(current_role.title, [])

        for guest_uc in candidates_for_this_role:
            # Prevent using the same unique ingredient object twice in the same arrangement
            # (Checks by unique ID)
            if any(g_info['obj'].id == guest_uc.id for g_info in current_arrangement_guests_info):
                logger.debug(f"Skipping {guest_uc.name} (ID: {guest_uc.id}) for role {current_role.title}: already in current arrangement.")
                continue

            # Skip if this guest_uc has no valid conformations (e.g., due to filtering issues in main.py)
            if not guest_uc.conformations or not os.path.exists(guest_uc.conformations):
                logger.debug(f"Skipping {guest_uc.name} (ID: {guest_uc.id}) for role {current_role.title}: no valid conformations file found at {guest_uc.conformations}.")
                continue

            # Read all valid conformations for this guest from its specific filtered XYZ file
            # These are expected to be sorted by energy (lower is better)
            all_conformations_data = read_xyz(guest_uc.conformations, logger) 

            for conf_idx, (total_atoms, conf_comment, conf_coords, conf_atom_types) in enumerate(all_conformations_data):
                # Parse energy from the comment line
                try:
                    conf_energy = float(conf_comment)
                except:
                    conf_energy = 0.0
                    logger.warning(f"Could not read docking score from conformer comment line for {guest_uc.conformations}")

                # Extract only guest coordinates and atom types from the combined host+guest conformation
                # This slicing assumes conf_coords contains host atoms first, then guest atoms.
                # If guest_uc.conformations points to guest-only file, this will produce empty array.
                guest_coords_only = conf_coords[host_atom_count:]
                guest_atom_types_only = conf_atom_types[host_atom_count:]

                # Skip if guest_coords_only is empty (likely means guest_uc.conformations was a guest-only file)
                if guest_coords_only.size == 0:
                    logger.debug(f"Skipping {guest_uc.name} (conf {conf_idx}) for role {current_role.title}: Guest coordinates are empty after slicing. Check if {guest_uc.conformations} contains host+guest.")
                    continue

                # Check for spatial overlap with already arranged guests
                overlap = False
                for existing_guest_info in current_arrangement_guests_info:
                    existing_guest_coords = existing_guest_info['coords']
                    existing_guest_atom_types = existing_guest_info['atom_types']

                    # Compare all atoms of the current guest with all atoms of each existing guest
                    for i, g1_coord in enumerate(guest_coords_only):
                        g1_atom_type = guest_atom_types_only[i]
                        for j, g2_coord in enumerate(existing_guest_coords):
                            g2_atom_type = existing_guest_atom_types[j]
                            distance = calculate_distance(g1_coord, g2_coord)
                            
                            # Overlap criterion: distance less than a fraction of sum of Van der Waals radii
                            vdw_sum = get_vdw_radius(g1_atom_type) + get_vdw_radius(g2_atom_type)
                            if distance < vdw_sum * 0.7:
                                overlap = True
                                # logger.debug(f"Overlap detected: {guest_uc.name} (conf {conf_idx}) and {existing_guest_info['obj'].name}. Dist: {distance:.2f} vs VdW sum: {vdw_sum:.2f} (threshold: {vdw_sum*0.7:.2f})")
                                break
                        if overlap:
                            break
                    if overlap:
                        break

                if not overlap:
                    # If no overlap, recursively try adding this guest to the arrangement
                    new_guest_info = {
                        'obj': guest_uc,
                        'coords': guest_coords_only,
                        'atom_types': guest_atom_types_only,
                        'conf_idx': conf_idx,
                        'energy': conf_energy
                    }
                    new_arrangement_guests_info = current_arrangement_guests_info + [new_guest_info]
                    
                    _recursively_arrange(
                        current_role_index + 1,
                        new_arrangement_guests_info,
                        satisfied_roles_count + 1,
                        current_energy_sum + conf_energy,
                    )
                    # If a valid, non-overlapping conformation is found for this guest,
                    # and recursive call finishes (doesn't raise MaxArrangementsReached),
                    # we can stop exploring other (worse energy) conformers for this specific guest.
                    break 

        # Option 2: Do NOT satisfy the current role (skip it)
        _recursively_arrange(
            current_role_index + 1,
            current_arrangement_guests_info, # Pass the same guests as before
            satisfied_roles_count,           # Do not increment satisfied_roles_count
            current_energy_sum               # Energy sum remains same
        )

    # Initial call to the recursive function
    try:
        # Start recursion with initial empty arrangement and 0 energy sum
        _recursively_arrange(0, [], 0, 0.0) 
    except MaxArrangementsReached:
        logger.info("Search terminated early as MAX_FINAL_ARRANGEMENTS limit was reached.")

    logger.info(f"All possibilities for unique arrangements satisfying given constraints have been exhausted.")
    # Format the return value based on the best_arrangement_info_per_signature
    return_list = []
    # Sort the final list by total_energy_sum (most negative first) for consistent output
    sorted_final_arrangements = sorted(best_arrangement_info_per_signature.values(), 
                                       key=lambda x: x['total_energy_sum'])

    for arr_item in sorted_final_arrangements:
        return_list.append({
            'guests': [
                {'name': g_info['obj'].name, 'role_title': g_info['obj'].role_title, 'conformation_idx': g_info['conf_idx']}
                for g_info in arr_item['guests_info']
            ],
            'satisfied_roles_count': arr_item['satisfied_roles_count']
        })

    if not return_list:
        logger.info("No valid arrangements could be found satisfying any roles within the given limits.")

    return return_list