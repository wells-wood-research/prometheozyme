import os
import numpy as np

from utils import get_vdw_radius, read_xyz, write_xyz, calculate_distance

# Constants for filtering and limiting output
MAX_FINAL_ARRANGEMENTS = 10  # Maximum number of unique arrangements to keep
RMSD_THRESHOLD = 2        # Angstroms - minimum RMSD for two arrangements to be considered distinct

# Custom exception for early termination
class MaxArrangementsReached(Exception):
    """Custom exception to signal that the maximum number of arrangements has been found."""
    pass

# Custom exception for early termination
class MaxArrangementsReached(Exception):
    """Custom exception to signal that the maximum number of arrangements has been found."""
    pass

def arrange_guests(roles, unique_guests_constraints, host_path, outdir, logger=None):
    """
    Finds a limited number of optimal and diverse arrangements of guests around a host,
    satisfying as many roles as possible according to their priority, and ensuring no
    spatial overlap between guests. Redundancy is reduced using RMSD filtering and
    by ensuring unique combinations of ingredients satisfying roles.
    The search stops when MAX_FINAL_ARRANGEMENTS are found.

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

    # Sort roles by priority (lowest priority value is most important)
    sorted_roles = sorted(roles, key=lambda r: r.priority)

    # This buffer will store the final optimal and diverse arrangements
    final_optimal_arrangements_buffer = []
    
    # This set will store unique signatures of arrangements (guest ID, role_title tuples)
    # to avoid permutations of the same ingredients satisfying the same roles.
    seen_arrangements_signatures = set() 

    max_satisfied_roles = -1

    # Read host coordinates and atom types once for combining later
    # Assuming read_xyz takes a logger and returns a list of (atom_count, comment, coords_np_array, atom_types_list)
    host_data = read_xyz(host_path, logger)
    if not host_data:
        logger.error(f"Could not read host file: {host_path}. Exiting arrangement process.")
        return []
    
    host_atom_count, host_comment, host_coords, host_atom_types = host_data[0]

    # Recursive helper function to explore arrangements
    def _recursively_arrange(
        current_role_index,
        current_arrangement_guests_info, # list of dicts: {'obj': guest_obj, 'coords': np.array, 'atom_types': list, 'conf_idx': int}
        satisfied_roles_count,
    ):
        nonlocal max_satisfied_roles

        # Base Case: All roles have been considered
        if current_role_index == len(sorted_roles):
            # If we found an arrangement that satisfies more roles, clear previous ones
            if satisfied_roles_count > max_satisfied_roles:
                max_satisfied_roles = satisfied_roles_count
                final_optimal_arrangements_buffer.clear() 
                seen_arrangements_signatures.clear() # Also clear seen signatures for higher optimality

            if satisfied_roles_count == max_satisfied_roles:
                # 1. Create the signature for the current arrangement
                # The signature identifies the unique set of guests and the roles they satisfied.
                current_signature_parts = []
                for g_info in current_arrangement_guests_info:
                    current_signature_parts.append((g_info['obj'].id, g_info['obj'].role_title))
                
                # Sort the parts to ensure a canonical representation regardless of order of addition
                current_signature = tuple(sorted(current_signature_parts)) 

                # Check if this combination of ingredients/roles has already been saved
                if current_signature in seen_arrangements_signatures:
                    logger.debug(f"Arrangement (roles: {satisfied_roles_count}) with identical ingredient/role combination already found. Skipping.")
                    return # Skip if identical signature found

                # 2. If signature is unique, then proceed with RMSD check for structural diversity
                
                # Combine all guest coordinates for the current arrangement
                if current_arrangement_guests_info:
                    new_arrangement_guest_coords_combined = np.concatenate([g_info['coords'] for g_info in current_arrangement_guests_info])
                else:
                    new_arrangement_guest_coords_combined = np.array([]) # No guests, empty arrangement

                is_structurally_unique = True
                
                # Only perform RMSD check if there are actual guest atoms to compare
                if new_arrangement_guest_coords_combined.size > 0:
                    for existing_arr_item in final_optimal_arrangements_buffer:
                        existing_arr_guests_info = existing_arr_item['guests_info']
                        
                        if existing_arr_guests_info:
                            existing_arr_guest_coords_combined = np.concatenate([g_info['coords'] for g_info in existing_arr_guests_info])
                        else:
                            existing_arr_guest_coords_combined = np.array([])

                        # Only compute RMSD if both arrangements have guests AND the same number of atoms
                        if (new_arrangement_guest_coords_combined.shape == existing_arr_guest_coords_combined.shape and
                            existing_arr_guest_coords_combined.size > 0): 
                            
                            # Simple RMSD calculation without alignment. Assumes corresponding atoms are in order.
                            rmsd = np.sqrt(np.sum((new_arrangement_guest_coords_combined - existing_arr_guest_coords_combined)**2) / new_arrangement_guest_coords_combined.shape[0])
                            
                            if rmsd < RMSD_THRESHOLD:
                                is_structurally_unique = False
                                logger.debug(f"Arrangement (roles: {satisfied_roles_count}) found structurally redundant (RMSD: {rmsd:.3f} < {RMSD_THRESHOLD}). Skipping.")
                                break
                    
                # Handle cases where both are empty arrangements (no guests)
                elif new_arrangement_guest_coords_combined.size == 0 and any(item['guests_info'].size == 0 for item in final_optimal_arrangements_buffer):
                    is_structurally_unique = False
                    logger.debug("Both current and an existing arrangement are empty, considered structurally redundant.")


                if is_structurally_unique:
                    if len(final_optimal_arrangements_buffer) < MAX_FINAL_ARRANGEMENTS:
                        final_optimal_arrangements_buffer.append({
                            'guests_info': current_arrangement_guests_info,
                            'satisfied_roles_count': satisfied_roles_count
                        })
                        seen_arrangements_signatures.add(current_signature) # Add signature to the set
                        logger.debug(f"Added unique arrangement (roles: {satisfied_roles_count}). Buffer size: {len(final_optimal_arrangements_buffer)}/{MAX_FINAL_ARRANGEMENTS}")
                        
                        # Check if buffer is now full, if so, raise exception to stop search
                        if len(final_optimal_arrangements_buffer) == MAX_FINAL_ARRANGEMENTS:
                            logger.info(f"Reached MAX_FINAL_ARRANGEMENTS ({MAX_FINAL_ARRANGEMENTS}). Stopping search early.")
                            raise MaxArrangementsReached
                    else:
                        logger.debug(f"Buffer is full ({MAX_FINAL_ARRANGEMENTS}). Skipping new unique arrangement (roles: {satisfied_roles_count}).")
            return # End of base case processing

        current_role = sorted_roles[current_role_index]
        
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
            # This file should contain the host + guest merged structures if main.py is updated
            all_conformations_data = read_xyz(guest_uc.conformations, logger) 

            for conf_idx, (total_atoms, conf_comment, conf_coords, conf_atom_types) in enumerate(all_conformations_data):
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
                                logger.debug(f"Overlap detected: {guest_uc.name} (conf {conf_idx}) and {existing_guest_info['obj'].name}. Dist: {distance:.2f} vs VdW sum: {vdw_sum:.2f} (threshold: {vdw_sum*0.7:.2f})")
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
                        'conf_idx': conf_idx
                    }
                    new_arrangement_guests_info = current_arrangement_guests_info + [new_guest_info]
                    _recursively_arrange(
                        current_role_index + 1,
                        new_arrangement_guests_info,
                        satisfied_roles_count + 1,
                    )

        # Option 2: Do NOT satisfy the current role (skip it)
        _recursively_arrange(
            current_role_index + 1,
            current_arrangement_guests_info, # Pass the same guests as before
            satisfied_roles_count,           # Do not increment satisfied_roles_count
        )

    # Initial call to the recursive function
    try:
        _recursively_arrange(0, [], 0)
    except MaxArrangementsReached:
        logger.info("Search terminated early as MAX_FINAL_ARRANGEMENTS limit was reached.")

    # Post-processing: Save the arrangements from the limited buffer
    if final_optimal_arrangements_buffer:
        logger.info(f"Saving {len(final_optimal_arrangements_buffer)} unique arrangement(s) satisfying {max_satisfied_roles} roles.")
        
        for i, arrangement_item in enumerate(final_optimal_arrangements_buffer):
            arrangement_guests_info = arrangement_item['guests_info']
            satisfied_roles_count = arrangement_item['satisfied_roles_count']

            combined_coords = host_coords.tolist()
            combined_atom_types = host_atom_types[:]
            total_atoms_in_arrangement = host_atom_count # Initialize with host's atom count
            
            # Construct a descriptive comment for the XYZ file
            arrangement_comment_parts = [f"Arrangement {i+1} ({satisfied_roles_count} roles satisfied)"]
            
            for guest_info in arrangement_guests_info:
                # Here, guest_info['obj'] is already the full unique_guests_constraints object
                # and guest_info['coords'] and guest_info['atom_types'] are the already extracted guest parts.
                
                combined_coords.extend(guest_info['coords'].tolist())
                combined_atom_types.extend(guest_info['atom_types'])
                total_atoms_in_arrangement += len(guest_info['coords'])
                arrangement_comment_parts.append(f"{guest_info['obj'].name}({guest_info['obj'].role_title}@{guest_info['conf_idx']})")

            arrangement_output_path = os.path.join(outdir, f"arrangement_{i+1}.xyz")
            # Use the imported write_xyz function from utils.py
            write_xyz(arrangement_output_path, " ".join(arrangement_comment_parts), np.array(combined_coords), combined_atom_types)
            logger.info(f"Arrangement {i+1} saved to {arrangement_output_path}")
    else:
        logger.info("No valid arrangements could be found satisfying any roles within the given limits.")

    # Format the return value to match the previous expected structure
    return_list = []
    for arr_buf_item in final_optimal_arrangements_buffer:
        return_list.append({
            'guests': [
                {'name': g_info['obj'].name, 'role_title': g_info['obj'].role_title, 'conformation_idx': g_info['conf_idx']}
                for g_info in arr_buf_item['guests_info']
            ],
            'satisfied_roles_count': arr_buf_item['satisfied_roles_count']
        })

    return return_list