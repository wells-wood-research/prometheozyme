from drOrca import make_orca_input
import os
import numpy as np
from utils import get_atom_count, read_xyz, calculate_center_of_mass, write_xyz #, add_dummy_atom # Placeholder for now

def calculate_charge(charges):
    return int(sum(charges))
    
def calculate_multiplicity(multiplicities):
    # Assume molecules are weakly or non-interacting
    spin = int(sum(m - 1 for m in multiplicities) / 2)
    multiplicity = 2*spin + 1
    return multiplicity

def optimise(arr, host_atom_count, ingredient_map, logger):
    original_xyz_path = arr["path"]
    arrName = os.path.splitext(os.path.basename(original_xyz_path))[0]
    orcaInputDir = os.path.join(os.path.dirname(original_xyz_path), arrName)
    os.makedirs(orcaInputDir, exist_ok=True)

    title_prefix = f"Optimisation of {arrName}"
    qmMethod = "XTB2"
    method = "Opt"
    inputFormat = "xyzfile"
    moleculeInfo = {"charge": calculate_charge([guest["obj"].charge for guest in arr["guests_info"]]),
                    "multiplicity": calculate_multiplicity([guest["obj"].multiplicity for guest in arr["guests_info"]])}
    parallelize = 8
    maxcore = 2500

    # Read the initial XYZ file
    xyz_data = read_xyz(original_xyz_path, logger)
    if not xyz_data:
        logger.error(f"Could not read XYZ data from {original_xyz_path}")
        return None # Or raise an error
    
    # Assuming single structure in the XYZ file for now
    original_atom_count, original_comment, initial_coordinates, initial_atom_types = xyz_data[0]

    orca_input_files = []
    
    # Placeholder for guest and host indices for COM calculation if needed
    # These would typically come from the input `arr` or be derived.
    # For now, let's assume guest_atom_indices and host_atom_indices cover all atoms of guest/host.
    # This part needs to be adapted based on how guest/host atoms are identified.
    # guest_atom_indices_relative_to_guest = list(range(guest_atom_count))
    # host_atom_indices_absolute = list(range(host_atom_count))


    # The main logic will iterate through constraints and generate files.
    # This simplified loop assumes one primary constraint drives the file generation.
    # In reality, you might have multiple guests or complex constraint interactions.
    
    # This is a simplified representation. The actual loop will be more complex
    # and depend on how constraints are structured in `arr["guests_info"]`.
    # We will iterate over each constraint to determine the generation strategy.

    for guest_info_idx, guest_info in enumerate(arr["guests_info"]):
        guest = guest_info["obj"]
        current_guest_atom_count = get_atom_count(ingredient_map[guest_info["obj"].name].path, logger)
        # This is the start index of the *first* guest. If multiple guests, this needs adjustment.
        # For now, assuming guest atoms are appended after host atoms.
        # guest_start_abs_idx = host_atom_count 

        # We need to correctly identify the range of atoms belonging to the current guest
        # and the host. This is crucial for COM calculations and indexing.
        # Let's assume guest_indices_for_com and host_indices_for_com are defined correctly
        # based on host_atom_count and accumulated guest atom counts.
        
        # This calculation of guest_start_abs_idx needs to be accurate for multiple guests.
        accumulated_guest_atom_count_so_far = sum(
            get_atom_count(ingredient_map[arr["guests_info"][i]["obj"].name].path, logger) 
            for i in range(guest_info_idx)
        )
        guest_start_abs_idx = host_atom_count + accumulated_guest_atom_count_so_far
        
        guest_atom_indices_for_com = list(range(guest_start_abs_idx, guest_start_abs_idx + current_guest_atom_count))
        host_atom_indices_for_com = list(range(host_atom_count))


        for constraint_idx, constraint in enumerate(guest.constraints):
            guestIdx_orig, guestType, hostIdx_orig, hostType, val = constraint

            # Work with copies for each constraint to handle potential DA additions independently
            current_xyz_coords = initial_coordinates.copy()
            current_atom_types = list(initial_atom_types) # Mutable copy
            
            # This will be the path to the XYZ file used by ORCA (might be a temporary one with DAs)
            processed_xyz_path = original_xyz_path 
            # Track DA additions for correct indexing and file writing
            guest_da_added = False
            host_da_added = False

            final_guest_indices = []
            final_host_indices = []

            # --- Guest Handling ---
            if guestType == 'com':
                guest_com = calculate_center_of_mass(current_xyz_coords, guest_atom_indices_for_com, current_atom_types)
                # Placeholder for add_dummy_atom call
                # temp_coords, temp_types, guest_DA_idx = add_dummy_atom(current_xyz_coords, current_atom_types, guest_com, "XG") 
                # current_xyz_coords, current_atom_types = temp_coords, temp_types
                # guest_da_added = True
                # final_guest_indices = [guest_DA_idx]
                
                # SIMULATED add_dummy_atom for guest
                current_xyz_coords = np.vstack([current_xyz_coords, guest_com])
                current_atom_types.append("XG") # XG for eXtra Guest (dummy)
                guest_DA_idx = len(current_atom_types) - 1 # 0-indexed
                final_guest_indices = [guest_DA_idx]
                guest_da_added = True
                logger.info(f"Guest COM at {guest_com}. Added dummy atom XG at index {guest_DA_idx}.")

            elif guestType == 'iter':
                # Adjust guestIdx_orig to be absolute indices. These are 0-indexed.
                final_guest_indices = [guest_start_abs_idx + idx for idx in guestIdx_orig]
            else:
                logger.error(f"Unsupported guestType: {guestType} for constraint {constraint_idx}")
                continue

            # --- Host Handling ---
            # Note: Host COM calculation uses original host indices, even if guest DA was added.
            if hostType == 'com':
                host_com = calculate_center_of_mass(current_xyz_coords, host_atom_indices_for_com, current_atom_types) # Pass current_atom_types which might include XG
                # Placeholder for add_dummy_atom call
                # temp_coords, temp_types, host_DA_idx = add_dummy_atom(current_xyz_coords, current_atom_types, host_com, "XH")
                # current_xyz_coords, current_atom_types = temp_coords, temp_types
                # host_da_added = True
                # final_host_indices = [host_DA_idx]

                # SIMULATED add_dummy_atom for host
                current_xyz_coords = np.vstack([current_xyz_coords, host_com])
                current_atom_types.append("XH") # XH for eXtra Host (dummy)
                host_DA_idx = len(current_atom_types) - 1 # 0-indexed
                final_host_indices = [host_DA_idx]
                host_da_added = True
                logger.info(f"Host COM at {host_com}. Added dummy atom XH at index {host_DA_idx}.")

            elif hostType == 'iter':
                # hostIdx_orig are already absolute 0-indexed indices.
                final_host_indices = list(hostIdx_orig)
            else:
                logger.error(f"Unsupported hostType: {hostType} for constraint {constraint_idx}")
                continue
            
            # --- Write temporary XYZ if DAs were added ---
            if guest_da_added or host_da_added:
                temp_xyz_filename = f"temp_{arrName}_g{guestType}_h{hostType}_c{constraint_idx}.xyz"
                processed_xyz_path = os.path.join(orcaInputDir, temp_xyz_filename)
                comment_for_temp_xyz = f"Temp XYZ for {arrName} with DAs (guest: {guestType}, host: {hostType}, constraint: {constraint_idx})"
                write_xyz(processed_xyz_path, comment_for_temp_xyz, current_xyz_coords, current_atom_types)
                logger.info(f"Written temporary XYZ with DAs to: {processed_xyz_path}")


            # --- File Generation Logic ---
            # Indices in final_guest_indices and final_host_indices are ALREADY absolute (0-indexed) 
            # and correct for the (potentially modified) structure.
            # ORCA requires 1-based indexing for atom specifications in constraint blocks.
            if guestType == 'iter' and hostType == 'iter':
                for g_idx_0based in final_guest_indices: 
                    for h_idx_0based in final_host_indices:
                        # Convert to 1-based for ORCA geom block and filenames/titles
                        g_idx_1based, h_idx_1based = g_idx_0based + 1, h_idx_0based + 1
                        geom = {"keep": [{"atoms": [g_idx_1based, h_idx_1based], "val": val}]}
                        orca_filename = f"opt_g{g_idx_1based}_h{h_idx_1based}_c{constraint_idx}.inp"
                        orcaInputPath = os.path.join(orcaInputDir, orca_filename)
                        title = f"{title_prefix} (constraint {constraint_idx}: guest_atom {g_idx_1based}, host_atom {h_idx_1based})"
                        make_orca_input(orcaInput=orcaInputPath, title=title, qmMethod=qmMethod, method=method,
                                        inputFormat=inputFormat, inputFile=processed_xyz_path, 
                                        moleculeInfo=moleculeInfo, parallelize=parallelize, maxcore=maxcore,
                                        geom=geom, qmmm=None, neb=None, scf=None, docker=None)
                        orca_input_files.append(orcaInputPath)
            elif guestType == 'iter' and hostType == 'com': # Host is DA
                h_da_idx_0based = final_host_indices[0] 
                for g_idx_0based in final_guest_indices: 
                    g_idx_1based, h_da_idx_1based = g_idx_0based + 1, h_da_idx_0based + 1
                    geom = {"keep": [{"atoms": [g_idx_1based, h_da_idx_1based], "val": val}]}
                    orca_filename = f"opt_g{g_idx_1based}_hDA{h_da_idx_1based}_c{constraint_idx}.inp"
                    orcaInputPath = os.path.join(orcaInputDir, orca_filename)
                    title = f"{title_prefix} (constraint {constraint_idx}: guest_atom {g_idx_1based}, host_DA {h_da_idx_1based})"
                    make_orca_input(orcaInput=orcaInputPath, title=title, qmMethod=qmMethod, method=method,
                                    inputFormat=inputFormat, inputFile=processed_xyz_path, 
                                    moleculeInfo=moleculeInfo, parallelize=parallelize, maxcore=maxcore,
                                    geom=geom, qmmm=None, neb=None, scf=None, docker=None)
                    orca_input_files.append(orcaInputPath)
            elif guestType == 'com' and hostType == 'iter': # Guest is DA
                g_da_idx_0based = final_guest_indices[0] 
                for h_idx_0based in final_host_indices: 
                    g_da_idx_1based, h_idx_1based = g_da_idx_0based + 1, h_idx_0based + 1
                    geom = {"keep": [{"atoms": [g_da_idx_1based, h_idx_1based], "val": val}]}
                    orca_filename = f"opt_gDA{g_da_idx_1based}_h{h_idx_1based}_c{constraint_idx}.inp"
                    orcaInputPath = os.path.join(orcaInputDir, orca_filename)
                    title = f"{title_prefix} (constraint {constraint_idx}: guest_DA {g_da_idx_1based}, host_atom {h_idx_1based})"
                    make_orca_input(orcaInput=orcaInputPath, title=title, qmMethod=qmMethod, method=method,
                                    inputFormat=inputFormat, inputFile=processed_xyz_path, 
                                    moleculeInfo=moleculeInfo, parallelize=parallelize, maxcore=maxcore,
                                    geom=geom, qmmm=None, neb=None, scf=None, docker=None)
                    orca_input_files.append(orcaInputPath)
            elif guestType == 'com' and hostType == 'com': # Both are DA
                g_da_idx_0based = final_guest_indices[0] 
                h_da_idx_0based = final_host_indices[0] 
                g_da_idx_1based, h_da_idx_1based = g_da_idx_0based + 1, h_da_idx_0based + 1
                geom = {"keep": [{"atoms": [g_da_idx_1based, h_da_idx_1based], "val": val}]}
                orca_filename = f"opt_gDA{g_da_idx_1based}_hDA{h_da_idx_1based}_c{constraint_idx}.inp"
                orcaInputPath = os.path.join(orcaInputDir, orca_filename)
                title = f"{title_prefix} (constraint {constraint_idx}: guest_DA {g_da_idx_1based}, host_DA {h_da_idx_1based})"
                make_orca_input(orcaInput=orcaInputPath, title=title, qmMethod=qmMethod, method=method,
                                inputFormat=inputFormat, inputFile=processed_xyz_path, 
                                moleculeInfo=moleculeInfo, parallelize=parallelize, maxcore=maxcore,
                                geom=geom, qmmm=None, neb=None, scf=None, docker=None)
                orca_input_files.append(orcaInputPath)

    if not orca_input_files:
        logger.warning(f"No ORCA input files were generated for {arrName} based on guest constraints.")
        return None # Or an empty list, depending on desired behavior
        
    # If only one file was generated and no dummy atoms were involved, 
    # could return the path to the .xyz, otherwise list of .inp
    # For now, always returning list of input files.
    return orca_input_files

# Note: The input_path argument for pull_backbone seems unused in the original code.
# I'll keep it in the signature if it's part of an interface but won't use it unless its purpose becomes clear.
def pull_backbone(arr, host_atom_count, ingredient_map, logger, input_path=None): # Added default for input_path
    original_xyz_path = arr["path"] # Assuming arr["path"] is the relevant XYZ for pulling
    arrName = os.path.splitext(os.path.basename(original_xyz_path))[0]
    # Create a subdirectory for pull outputs, possibly different from "opt" outputs
    pullOrcaInputDir = os.path.join(os.path.dirname(original_xyz_path), f"{arrName}_pull") 
    os.makedirs(pullOrcaInputDir, exist_ok=True)

    title_prefix = f"Pulling backbone for {arrName}"
    qmMethod = "XTB2" # Assuming same method as optimise
    method = "Opt"    # Assuming same method as optimise
    inputFormat = "xyzfile"
    moleculeInfo = {"charge": calculate_charge([guest["obj"].charge for guest in arr["guests_info"]]),
                    "multiplicity": calculate_multiplicity([guest["obj"].multiplicity for guest in arr["guests_info"]])}
    parallelize = 8   # Assuming same settings
    maxcore = 2500    # Assuming same settings

    xyz_data = read_xyz(original_xyz_path, logger)
    if not xyz_data:
        logger.error(f"Could not read XYZ data from {original_xyz_path} for pull_backbone")
        return None
    
    original_atom_count, original_comment, initial_coordinates, initial_atom_types = xyz_data[0]
    orca_input_files = []

    for guest_info_idx, guest_info in enumerate(arr["guests_info"]):
        guest = guest_info["obj"]
        current_guest_atom_count = get_atom_count(ingredient_map[guest_info["obj"].name].path, logger)
        
        accumulated_guest_atom_count_so_far = sum(
            get_atom_count(ingredient_map[arr["guests_info"][i]["obj"].name].path, logger) 
            for i in range(guest_info_idx)
        )
        guest_start_abs_idx = host_atom_count + accumulated_guest_atom_count_so_far
        
        guest_atom_indices_for_com = list(range(guest_start_abs_idx, guest_start_abs_idx + current_guest_atom_count))
        host_atom_indices_for_com = list(range(host_atom_count))

        for constraint_idx, constraint in enumerate(guest.constraints):
            guestIdx_orig, guestType, hostIdx_orig, hostType, val_keep = constraint # val_keep is for the "keep" constraint
            val_pots = 0.5 # Potentials have a constant force, as per original logic

            current_xyz_coords = initial_coordinates.copy()
            current_atom_types = list(initial_atom_types)
            processed_xyz_path = original_xyz_path
            guest_da_added = False
            host_da_added = False
            final_guest_indices = []
            final_host_indices = []

            # --- Guest Handling (same as in optimise) ---
            if guestType == 'com':
                guest_com = calculate_center_of_mass(current_xyz_coords, guest_atom_indices_for_com, current_atom_types)
                current_xyz_coords = np.vstack([current_xyz_coords, guest_com])
                current_atom_types.append("XG")
                guest_DA_idx = len(current_atom_types) - 1
                final_guest_indices = [guest_DA_idx]
                guest_da_added = True
                logger.info(f"Pull: Guest COM at {guest_com}. Added XG at index {guest_DA_idx}.")
            elif guestType == 'iter':
                final_guest_indices = [guest_start_abs_idx + idx for idx in guestIdx_orig]
            else:
                logger.error(f"Pull: Unsupported guestType: {guestType} for constraint {constraint_idx}")
                continue

            # --- Host Handling (same as in optimise) ---
            if hostType == 'com':
                host_com = calculate_center_of_mass(current_xyz_coords, host_atom_indices_for_com, current_atom_types)
                current_xyz_coords = np.vstack([current_xyz_coords, host_com])
                current_atom_types.append("XH")
                host_DA_idx = len(current_atom_types) - 1
                final_host_indices = [host_DA_idx]
                host_da_added = True
                logger.info(f"Pull: Host COM at {host_com}. Added XH at index {host_DA_idx}.")
            elif hostType == 'iter':
                final_host_indices = list(hostIdx_orig)
            else:
                logger.error(f"Pull: Unsupported hostType: {hostType} for constraint {constraint_idx}")
                continue
            
            if guest_da_added or host_da_added:
                temp_xyz_filename = f"temp_pull_{arrName}_g{guestType}_h{hostType}_c{constraint_idx}.xyz"
                processed_xyz_path = os.path.join(pullOrcaInputDir, temp_xyz_filename) # Use pullOrcaInputDir
                comment_for_temp_xyz = f"Temp XYZ for pull {arrName} with DAs (g: {guestType}, h: {hostType}, c: {constraint_idx})"
                write_xyz(processed_xyz_path, comment_for_temp_xyz, current_xyz_coords, current_atom_types)
                logger.info(f"Pull: Written temporary XYZ with DAs to: {processed_xyz_path}")

            # --- File Generation Logic (with "keep" and "pots") ---
            def generate_pull_files(g_indices_0based, h_indices_0based, file_suffix_g, file_suffix_h):
                for g_idx_0based in g_indices_0based:
                    for h_idx_0based in h_indices_0based:
                        g_idx_1based, h_idx_1based = g_idx_0based + 1, h_idx_0based + 1
                        
                        geom = {
                            "keep": [{"atoms": [g_idx_1based, h_idx_1based], "val": val_keep}],
                            "pots": [{"atoms": [g_idx_1based, h_idx_1based], "val": val_pots}] 
                        }
                        # Note: Original pull_backbone used guest_start_abs_idx for pots' g_idx.
                        # Here, we use the same g_idx for keep and pots, which seems more consistent
                        # with constraining a specific guest atom/DA relative to a host atom/DA.
                        # If the original intent for pots was *always* the first atom of the guest, 
                        # that logic would need to be re-introduced carefully.
                        # For now, pots apply to the same pair as keep.

                        orca_filename = f"pull_g{file_suffix_g(g_idx_1based)}_h{file_suffix_h(h_idx_1based)}_c{constraint_idx}.inp"
                        orcaInputPath = os.path.join(pullOrcaInputDir, orca_filename) # Use pullOrcaInputDir
                        title = f"{title_prefix} (c {constraint_idx}: g {file_suffix_g(g_idx_1based)}, h {file_suffix_h(h_idx_1based)})"
                        
                        make_orca_input(orcaInput=orcaInputPath, title=title, qmMethod=qmMethod, method=method,
                                        inputFormat=inputFormat, inputFile=processed_xyz_path, 
                                        moleculeInfo=moleculeInfo, parallelize=parallelize, maxcore=maxcore,
                                        geom=geom, qmmm=None, neb=None, scf=None, docker=None)
                        orca_input_files.append(orcaInputPath)

            if guestType == 'iter' and hostType == 'iter':
                generate_pull_files(final_guest_indices, final_host_indices, lambda g: str(g), lambda h: str(h))
            elif guestType == 'iter' and hostType == 'com': # Host is DA
                h_da_idx_0based = final_host_indices[0]
                generate_pull_files(final_guest_indices, [h_da_idx_0based], lambda g: str(g), lambda h: f"DA{h}")
            elif guestType == 'com' and hostType == 'iter': # Guest is DA
                g_da_idx_0based = final_guest_indices[0]
                generate_pull_files([g_da_idx_0based], final_host_indices, lambda g: f"DA{g}", lambda h: str(h))
            elif guestType == 'com' and hostType == 'com': # Both are DA
                g_da_idx_0based = final_guest_indices[0]
                h_da_idx_0based = final_host_indices[0]
                generate_pull_files([g_da_idx_0based], [h_da_idx_0based], lambda g: f"DA{g}", lambda h: f"DA{h}")

    if not orca_input_files:
        logger.warning(f"Pull: No ORCA input files were generated for {arrName} based on guest constraints.")
        return None 
        
    return orca_input_files