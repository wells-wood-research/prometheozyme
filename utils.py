import numpy as np

# Define atomic masses for center of mass calculation
atomic_masses = {
        'H': 1.00794, 'He': 4.002602, 'C': 12.0107, 'N': 14.0067, 'O': 15.9994,
        'F': 18.998403, 'P': 30.973762, 'S': 32.065, 'Cl': 35.453,
        # Add more elements as needed
    }

# Define Van der Waals radii (simplified for example - use more accurate values if available)
vdw_radii = {
    'C': 1.7, 'H': 1.2, 'N': 1.55, 'O': 1.52, 'S': 1.8
}

def get_mass(atom):
    """Return the atomic mass of an element."""
    return atomic_masses.get(atom, 12.01)  # Default to carbon if unknown

def get_vdw_radius(atom):
    """Return the Van der Waals radius of an element."""
    return vdw_radii.get(atom, 1.7)  # Default to carbon if unknown

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

def write_xyz(filename, comment, coords, atom_types):
    """Write coordinates to an XYZ file."""
    with open(filename, 'w') as f:
        f.write(f"{len(coords)}\n")
        f.write(f"{comment}\n")
        for atom, (x, y, z) in zip(atom_types, coords):
            f.write(f"{atom} {x:.6f} {y:.6f} {z:.6f}\n")

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