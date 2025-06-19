from drOrca import make_orca_input
import os
from utils import get_atom_count, read_xyz, calculate_center_of_mass, calculate_distance, add_dummy_atom_to_xyz, write_xyz

def calculate_charge(charges):
    return int(sum(charges))
    
def calculate_multiplicity(multiplicities):
    # Assume molecules are weakly or non-interacting
    spin = int(sum(m - 1 for m in multiplicities) / 2)
    multiplicity = 2*spin + 1
    return multiplicity

def find_indices_closest_to_val(conformation_coords, absolute_guest_indices, host_indices, val):
    if not absolute_guest_indices or not host_indices:
        return (None, None), None
    min_diff = float('inf')
    best_pair = None
    best_distance = None

    for g_idx in absolute_guest_indices:
        guest_coord = conformation_coords[g_idx]
        for h_idx in host_indices:
            host_coord = conformation_coords[h_idx]
            distance = calculate_distance(guest_coord, host_coord)
            diff = abs(distance - val)
            if diff < min_diff:
                min_diff = diff
                best_pair = (g_idx, h_idx)
                best_distance = distance

    return best_pair, best_distance if best_pair is not None else (None, None)

def process_constraint_com(conformation_coords, atom_types, absolute_guest_indices, host_indices):
    """Evaluate constraints using center of mass distance."""
    guest_com = calculate_center_of_mass(conformation_coords, [i for i in absolute_guest_indices], atom_types)
    host_com = calculate_center_of_mass(conformation_coords, host_indices, atom_types)
    
    if guest_com is None or host_com is None:
        return (None, None), conformation_coords, atom_types
    
    updated_xyz_coords, updated_atom_types, guest_da_index = add_dummy_atom_to_xyz(conformation_coords, atom_types, guest_com)
    updated_xyz_coords, updated_atom_types, host_da_index = add_dummy_atom_to_xyz(updated_xyz_coords, updated_atom_types, host_com)
    
    return (guest_da_index, host_da_index), updated_xyz_coords, updated_atom_types

def process_constraint_mixed(path, original_comment, conformation_coords, current_atom_types, absolute_guest_indices, guestType, host_indices, hostType, val):
    """Evaluate constraints for mixed iter and com types."""
    if guestType == "iter" and hostType == "com":
        host_com = calculate_center_of_mass(conformation_coords, host_indices, current_atom_types)
        if host_com is None:
            return (None, None), conformation_coords, current_atom_types
        updated_xyz_coords, updated_atom_types, host_da_index = add_dummy_atom_to_xyz(conformation_coords, current_atom_types, host_com)
        result = find_indices_closest_to_val(updated_xyz_coords, absolute_guest_indices, [host_da_index], val)
        if result[0] is None:
            return (None, None), conformation_coords, current_atom_types
        (g_idx, h_idx), best_distance = result
        return (g_idx, h_idx), updated_xyz_coords, updated_atom_types
    
    elif guestType == "com" and hostType == "iter":
        guest_com = calculate_center_of_mass(conformation_coords, absolute_guest_indices, current_atom_types)
        if guest_com is None:
            return (None, None), conformation_coords, current_atom_types
        updated_xyz_coords, updated_atom_types, guest_da_index = add_dummy_atom_to_xyz(conformation_coords, current_atom_types, guest_com)
        result = find_indices_closest_to_val(updated_xyz_coords, [guest_da_index], host_indices, val)
        if result[0] is None:
            return (None, None), conformation_coords, current_atom_types
        (g_idx, h_idx), best_distance = result
        return (g_idx, h_idx), updated_xyz_coords, updated_atom_types

def optimise(arr, host_atom_count, ingredient_map, backbone, logger):
    path = arr["path"]
    arrName = os.path.splitext(os.path.basename(path))[0]
    orcaInputDir = os.path.join(os.path.dirname(path), arrName)
    orcaInputDir_Opt = os.path.join(orcaInputDir, "opt")
    orcaInputDir_Pull = os.path.join(orcaInputDir, "pull")
    os.makedirs(orcaInputDir_Opt, exist_ok=True)
    os.makedirs(orcaInputDir_Pull, exist_ok=True)
    orcaInput_Opt = os.path.join(orcaInputDir_Opt, "opt.inp")
    orcaInput_Pull = os.path.join(orcaInputDir_Pull, "pull.inp")

    title = f"Initial optimisation of {arrName}:\n # {arr['desc']}"
    
    qmMethod = "XTB2"
    method = "Opt"
    inputFormat = "xyzfile"

    moleculeInfo = {"charge": calculate_charge([guest["obj"].charge for guest in arr["guests_info"]]),
                    "multiplicity": calculate_multiplicity([guest["obj"].multiplicity for guest in arr["guests_info"]])}

    parallelize = 8
    maxcore = 2500

    xyz_data = read_xyz(path, logger)
    if not xyz_data:
        logger.error(f"Could not read XYZ data from {path}")
        return None
    original_atom_count, original_comment, initial_coordinates, initial_atom_types = xyz_data[0]

    geom = {"keep": []}
    geom_pots = {"pots": []}
    
    current_xyz_coords = initial_coordinates.copy()
    current_atom_types = list(initial_atom_types)
    
    accumulated_guest_atom_count = 0
    for guest_info in arr["guests_info"]:
        guest = guest_info["obj"]
        current_guest_atom_count = get_atom_count(ingredient_map[guest_info["obj"].name].path)
        current_guest_absolute_start_index = host_atom_count + accumulated_guest_atom_count
        if backbone:
            backbone_idx = current_guest_absolute_start_index + 1
        else:
            backbone_idx = current_guest_absolute_start_index
        for constraint_idx, constraint in enumerate(guest.constraints):
            guestIdx_orig, guestType, hostIdx_orig, hostType, val = constraint
            absolute_guest_indices = [current_guest_absolute_start_index + g_idx for g_idx in guestIdx_orig]

            if guestType == "iter" and hostType == "iter":
                result = find_indices_closest_to_val(current_xyz_coords, absolute_guest_indices, hostIdx_orig, val)
                if result[0] is None:
                    logger.warning(f"Could not find suitable atom pair for constraint {constraint_idx} in guest {guest.name}")
                    continue
                (g_idx, h_idx), best_distance = result
                geom["keep"].append({"atoms": [g_idx, h_idx], "val": val})
                geom_pots["pots"].append({"atoms": [backbone_idx, h_idx], "val": 0.5})

            elif guestType == "com" and hostType == "com":
                result = process_constraint_com(current_xyz_coords, current_atom_types, absolute_guest_indices, hostIdx_orig)
                if result[0] == (None, None):
                    logger.warning(f"Could not process 'com' constraint {constraint_idx} in guest {guest.name}")
                    continue
                (g_idx, h_idx), current_xyz_coords, current_atom_types = result
                geom["keep"].append({"atoms": [g_idx, h_idx], "val": val})
                geom_pots["pots"].append({"atoms": [backbone_idx, h_idx], "val": 0.5})

            elif (guestType == "iter" and hostType == "com") or (guestType == "com" and hostType == "iter"):
                result = process_constraint_mixed(path, original_comment, current_xyz_coords, current_atom_types, absolute_guest_indices, guestType, hostIdx_orig, hostType, val)
                if result[0] == (None, None):
                    logger.warning(f"Could not process mixed constraint {constraint_idx} in guest {guest.name}")
                    continue
                (g_idx, h_idx), current_xyz_coords, current_atom_types = result
                geom["keep"].append({"atoms": [g_idx, h_idx], "val": val})
                geom_pots["pots"].append({"atoms": [backbone_idx, h_idx], "val": 0.5})

            else:
                logger.error(f"guestType in constraints for role {arr['desc']} cannot be {guestType}; only allowed options are 'iter' and 'com'!")
        accumulated_guest_atom_count += current_guest_atom_count

    # Write the final XYZ file with all dummy atoms
    updated_xyz_path = os.path.join(orcaInputDir_Opt, f"{arrName}_with_dummies.xyz")
    write_xyz(updated_xyz_path, f"{original_comment} (with Dummy Atoms)", current_xyz_coords, current_atom_types)

    # Use the updated XYZ file for optimization
    make_orca_input(orcaInput=orcaInput_Opt,
                    title=title,
                    qmMethod=qmMethod,
                    method=method,
                    inputFormat=inputFormat,
                    inputFile=updated_xyz_path,
                    moleculeInfo=moleculeInfo,
                    parallelize=parallelize,
                    maxcore=maxcore,
                    qmmm=None,
                    geom=geom,
                    neb=None,
                    scf=None,
                    docker=None)
    arr_optimised = orcaInput_Opt.replace(".inp", ".xyz")

    try:
        geom = geom | geom_pots
    except:
        try:
            geom = {**geom, **geom_pots}
        except:
            logger.error("Cannot merge geometry dictionaries for ORCA constraints")
    make_orca_input(orcaInput=orcaInput_Pull,
                    title=title,
                    qmMethod=qmMethod,
                    method=method,
                    inputFormat=inputFormat,
                    inputFile=arr_optimised,
                    moleculeInfo=moleculeInfo,
                    parallelize=parallelize,
                    maxcore=maxcore,
                    qmmm=None,
                    geom=geom,
                    neb=None,
                    scf=None,
                    docker=None)
    arr_pulled = orcaInput_Pull.replace(".inp", ".xyz")

    return arr_optimised, arr_pulled