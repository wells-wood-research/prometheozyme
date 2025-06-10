import numpy as np
import os

def read_xyz(file_path, logger):
    """Read an XYZ file and return a list of (atom_count, comment, coordinates, atom_types) where coordinates is a numpy array."""
    with open(file_path, 'r') as f:
        lines = [line.strip() for line in f.readlines()]
    
    structures = []
    i = 0
    while i < len(lines):
        # Skip empty lines
        while i < len(lines) and (not lines[i] or lines[i].isspace()):
            i += 1
        if i >= len(lines):
            break
        
        try:
            # Read atom count (must be a valid integer)
            atom_count = int(lines[i])
            i += 1
            if i >= len(lines):
                logger.warning(f"Unexpected end of file at line {i}: no comment line after atom count")
                break
            # Read comment line
            comment = lines[i]
            i += 1
            if i >= len(lines):
                logger.warning(f"Unexpected end of file at line {i}: no coords after comment")
                break
            
            # Read coordinates and atom types
            coords = []
            atom_types = []
            for j in range(i, min(i + atom_count, len(lines))):
                parts = lines[j].split()
                if len(parts) >= 4:
                    try:
                        x, y, z = map(float, parts[1:4])
                        coords.append([x, y, z])
                        atom_types.append(parts[0])
                    except ValueError:
                        logger.warning(f"Skipping invalid coordinate line at {j+1}")
                        continue
                else:
                    logger.warning(f"Skipping malformed line at {j+1}")
                    continue
            
            # Only append if we have valid coordinates
            if coords and len(coords) == atom_count:
                coordinates = np.array(coords)
                structures.append((atom_count, comment, coordinates, atom_types))
            else:
                logger.warning(f"Skipping structure at line {i-atom_count-1}: expected {atom_count} atoms, found {len(coords)}")
            
            i += atom_count
        except ValueError:
            logger.warning(f"Invalid atom count at line {i+1}, skipping to next line")
            i += 1
            continue
    
    return structures

def get_atom_count(path):
    """Read host XYZ file and return number of host atoms."""
    with open(path, 'r') as f:
        lines = f.readlines()
        if len(lines) < 1:
            return 0
        return int(lines[0].strip())  # First line is atom count

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

def filter_conformations(merged_path, host_path, id, name, role, constraints, logger):
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
    
    # Filter conformations
    valid_structures = []
    for atom_count, comment, coordinates, atom_types in structures:
        # Check if all constraints are satisfied
        all_constraints_satisfied = True
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
    if len(valid_structures) != 0:
        # Write filtered conformations to new file to ensure docked output that might be needed later is not affected
        filtered_path = os.path.join(os.path.dirname(merged_path), f"{name}_{role}_{id}.xyz")
        with open(filtered_path, 'w') as f:
            for atom_count, comment, coordinates, atom_types in valid_structures:
                f.write(f"{atom_count}\n")
                f.write(f"{comment}\n")
                for atom_type, coord in zip(atom_types, coordinates):
                    f.write(f"{atom_type} {coord[0]:27.17f} {coord[1]:27.17f} {coord[2]:27.17f}\n")
    
    return valid_structures, None if not filtered_path else filtered_path
