import os
import numpy as np

from utils import read_xyz, write_xyz

def get_atom_count(path, logger=None):
    """Read molecule file (XYZ or PDB) and return number of atoms."""
    _, ext = os.path.splitext(path)
    ext = ext.lower()

    if ext == '.xyz':
        try:
            with open(path, 'r') as f:
                lines = f.readlines()
                if len(lines) < 1:
                    if logger:
                        logger.warning(f"XYZ file {path} is empty.")
                    return 0
                return int(lines[0].strip())  # First line is atom count
        except FileNotFoundError:
            if logger:
                logger.error(f"XYZ file not found: {path}")
            return 0
        except ValueError:
            if logger:
                logger.error(f"Invalid atom count in XYZ file: {path}")
            return 0
        except Exception as e:
            if logger:
                logger.error(f"Error reading XYZ atom count from {path}: {e}")
            return 0
    elif ext == '.pdb':
        try:
            from pdbUtils import pdb2df
            df = pdb2df(path)
            return len(df) # Number of rows in DataFrame is atom count
        except FileNotFoundError:
            if logger:
                logger.error(f"PDB file not found: {path}")
            return 0
        except ImportError:
            if logger:
                logger.error("pdbUtils.pdb2df is required to read PDB files but could not be imported.")
            return 0
        except Exception as e:
            if logger:
                logger.error(f"Error reading PDB atom count from {path}: {e}")
            return 0
    else:
        if logger:
            logger.warning(f"Unsupported file format for atom count: {ext} for file {path}")
        return 0
    
def calculate_distance(coord1, coord2):
    """Calculate Euclidean distance between two 3D coordinates."""
    return np.sqrt(np.sum((coord1 - coord2) ** 2))

def calculate_center_of_mass(coordinates, indices, atom_types):
    """Calculate center of mass for given indices, weighted by atomic masses."""
    if not indices:
        return None
    
    # Atomic masses in atomic mass units (u) for common elements
    atomic_masses = {
        'H': 1.00794, 'He': 4.002602, 'C': 12.0107, 'N': 14.0067, 'O': 15.9994,
        'F': 18.998403, 'P': 30.973762, 'S': 32.065, 'Cl': 35.453,
        # Add more elements as needed
    }
    
    selected_coords = np.array([coordinates[i] for i in indices])
    masses = np.array([atomic_masses.get(atom_types[i], 1.0) for i in indices])  # Default to 1.0 if unknown
    
    # Weighted average: sum(mass * coord) / sum(mass)
    weighted_coords = selected_coords * masses[:, np.newaxis]
    center_of_mass = np.sum(weighted_coords, axis=0) / np.sum(masses)
    
    return center_of_mass

def evaluate_constraint_iter(conformation_coords, guest_indices, host_indices, val, host_atom_count):
    """Evaluate constraints iteratively (any guest-host atom pair satisfying keeps the conformation)."""
    for g_idx in guest_indices:
        guest_coord = conformation_coords[g_idx + host_atom_count]
        for h_idx in host_indices:
            host_coord = conformation_coords[h_idx]
            distance = calculate_distance(guest_coord, host_coord)
            # If any pair satisfies the constraint, return True
            if distance <= val + 0.5: # val - 0.5 <= distance <= val + 0.5
                return True
    return False

def evaluate_constraint_com(conformation_coords, guest_indices, host_indices, val, host_atom_count, atom_types):
    """Evaluate constraints using center of mass distance."""
    guest_com = calculate_center_of_mass(conformation_coords, [i + host_atom_count for i in guest_indices], atom_types)
    host_com = calculate_center_of_mass(conformation_coords, host_indices, atom_types)
    
    if guest_com is None or host_com is None:
        return False
    
    distance = calculate_distance(guest_com, host_com)
    return distance <= val + 0.5 # val - 0.5 <= distance <= val + 0.5

def evaluate_constraint_mixed(conformation_coords, guest_indices, guestType, host_indices, hostType, val, host_atom_count):
    """Evaluate constraints for mixed iter and com types."""
    if guestType == "iter" and hostType == "com":
        host_com = calculate_center_of_mass(conformation_coords, host_indices)
        if host_com is None:
            return False
        for g_idx in guest_indices:
            guest_coord = conformation_coords[g_idx + host_atom_count]
            distance = calculate_distance(guest_coord, host_com)
            if distance <= val + 0.5: # val - 0.5 <= distance <= val + 0.5
                return True
        return False
    elif guestType == "com" and hostType == "iter":
        guest_com = calculate_center_of_mass(conformation_coords, [i + host_atom_count for i in guest_indices])
        if guest_com is None:
            return False
        for h_idx in host_indices:
            host_coord = conformation_coords[h_idx]
            distance = calculate_distance(host_coord, guest_com)
            if distance <= val + 0.5: # val - 0.5 <= distance <= val + 0.5
                return True
        return False

def evaluate_constraint(coordinates, atom_types, guest_indices, guestType, host_indices, hostType, val, host_atom_count, logger):
    """Evaluate constraints based on guestType and hostType."""
    if guestType == "iter" and hostType == "iter":
        return evaluate_constraint_iter(coordinates, guest_indices, host_indices, val, host_atom_count)
    elif guestType == "com" and hostType == "com":
        return evaluate_constraint_com(coordinates, guest_indices, host_indices, val, host_atom_count, atom_types)
    elif (guestType == "iter" and hostType == "com") or (guestType == "com" and hostType == "iter"):
        return evaluate_constraint_mixed(coordinates, guest_indices, guestType, host_indices, hostType, val, host_atom_count, atom_types)
    else:
        logger.error(f"Invalid constraint types: guestType={guestType}, hostType={hostType}")
        return False

def evaluate_backbone_out(dir, coordinates, atom_types, host_atom_count):
    cb = np.array(coordinates[dir[0] + host_atom_count])
    cg = np.array(coordinates[dir[1] + host_atom_count])
    sidechain_vec = cg - cb

    host_coords = np.array([coord for atom, coord in zip(atom_types[:host_atom_count], coordinates[:host_atom_count]) if atom.upper() != 'H'])

    option = "nearest" # doesn't work very well with 'center'
    if option == "center":
        # Option 1: Use center of host
        host_center = np.mean(host_coords, axis=0)
        host_direction = host_center - cb
    elif option == "nearest":
        # Option 2: Use nearest host atom
        distances = np.linalg.norm(host_coords - cb, axis=1)
        nearest_host = host_coords[np.argmin(distances)]
        host_direction = nearest_host - cb
        
    dot_product = np.dot(sidechain_vec, host_direction)

    return dot_product < 0 # True if guest faces away

def filter_conformations(merged_path, host_path, name, role, indices, constraints, evaluate_backbone, logger=None):
    """Filter conformations in XYZ file based on multiple distance constraints."""
    # Get atom counts
    total_atoms = get_atom_count(merged_path)
    host_atom_count = get_atom_count(host_path)
    if host_atom_count == 0:
        logger.error(f"Invalid host file at {host_path}")
        return
    guest_atom_count = total_atoms - host_atom_count
    
    # Read all conformations
    structures = read_xyz(merged_path, logger)
    if not structures:
        logger.error(f"No conformations found in {merged_path}")
        return
    
    # Get atom indices for direction (backbone orientation) check
    dir = getattr(indices, "dir", [0,1])
    
    # Filter conformations
    valid_structures = []
    for atom_count, comment, coordinates, atom_types in structures:
        all_constraints_satisfied = True
        # Check if all constraints are satisfied
        if evaluate_backbone:
            if not evaluate_backbone_out(dir, coordinates, atom_types, host_atom_count):
                all_constraints_satisfied = False
                continue
        if constraints:
            for constraint in constraints:
                guest_indices, guestType, host_indices, hostType, val = constraint
                if guest_indices == "all":
                    guest_indices = list(range(0, guest_atom_count))
                if host_indices == "all":
                    host_indices = list(range(0, host_atom_count))
                if not evaluate_constraint(coordinates, atom_types, guest_indices, guestType, host_indices, hostType, val, host_atom_count, logger):
                    all_constraints_satisfied = False
                    break
        
        if all_constraints_satisfied:
            valid_structures.append((atom_count, comment, coordinates, atom_types))

    logger.info(f"Filtered {len(structures)} conformations to {len(valid_structures)} valid conformations")
    
    # TODO Adapt this to reuse write_xyz from utils
    filtered_path = None
    if len(valid_structures) != 0:
        # Write filtered conformations to new file to ensure docked output that might be needed later is not affected
        filtered_path = os.path.join(os.path.dirname(merged_path), f"{name}_{role}.xyz")
        with open(filtered_path, 'w') as f:
            for atom_count, comment, coordinates, atom_types in valid_structures:
                write_xyz(f, comment, coordinates, atom_types)
    
    return valid_structures, filtered_path
