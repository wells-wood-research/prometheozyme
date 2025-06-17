import numpy as np
import io
import os

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

def read_xyz(file_path, logger=None):
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

def merge_xyz(receptor_path, ligand_path, output_path):
    """Merge receptor XYZ with each structure in a multi-XYZ ligand file."""
    # Read receptor (single structure)
    receptor_structures = read_xyz(receptor_path)
    if not receptor_structures:
        raise ValueError("Receptor file is empty or invalid")
    receptor_count, receptor_comment, receptor_coords, receptor_atom_types = receptor_structures[0]
    
    # Read ligand (potentially multi-XYZ)
    ligand_structures = read_xyz(ligand_path)
    if not ligand_structures:
        raise ValueError("Ligand file is empty or invalid")
    
    # Write merged multi-XYZ file
    with open(output_path, 'w') as f:
        for ligand_count, ligand_comment, ligand_coords, ligand_atom_types in ligand_structures:
            # Create comment for merged structure
            merged_comment = f"Merged from {receptor_path} and {ligand_path} (structure)"
            merged_coordinates = np.concatenate([receptor_coords, ligand_coords], axis=0)
            merged_atom_types = receptor_atom_types + ligand_atom_types
            write_xyz(f, merged_comment, merged_coordinates, merged_atom_types)

def write_xyz(file, comment, coords, atom_types):
    """Write coordinates to an XYZ file."""
    def write_to_file(f):
        """Write the XYZ content to the given file object."""
        f.write(f"{len(coords)}\n")
        f.write(f"{comment}\n")
        for atom, (x, y, z) in zip(atom_types, coords):
            f.write(f"{atom} {x:27.17f} {y:27.17f} {z:27.17f}\n")
        f.write("\n")

    # Handle file object or path
    if isinstance(file, io.IOBase):
        write_to_file(file)
    elif isinstance(file, str):
        with open(file, 'w') as f:
            write_to_file(f)

# Modified get_atom_count to handle both XYZ and PDB
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

def read_score(filepath, logger=None):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    # Find the start of the results table
    start_index = None
    for i, line in enumerate(lines):
        if line.strip().startswith("mode |  affinity"):
            start_index = i + 3  # data starts 3 lines after the header line
            break

    # Parse the affinity column
    affinities = []
    if start_index is not None:
        for line in lines[start_index:]:
            parts = line.split()
            if len(parts) >= 2:
                try:
                    affinity = float(parts[1])  # Second column is "affinity"
                    affinities.append(affinity)
                except ValueError:
                    continue  # Skip lines that don't contain floats in expected place
    else:
        logger.error(f"No scores found in {filepath}!")

    # Convert to NumPy array or DataFrame
    affinity_array = np.array(affinities)

    return affinity_array

def append_scores(xyz_file, scores_file, logger=None):
    structures = read_xyz(xyz_file, logger=logger)
    scores = read_score(scores_file, logger=logger)

    if len(structures) != len(scores):
        if logger:
            logger.error(f"Mismatch: {len(structures)} structures vs {len(scores)} scores in {xyz_file}")
        else:
            raise ValueError(f"Mismatch: {len(structures)} structures vs {len(scores)} scores")

    temp_output = xyz_file + ".tmp"

    with open(temp_output, 'w') as f:
        for i, (atom_count, comment, coordinates, atom_types) in enumerate(structures):
            score = scores[i]
            new_comment = f"{score:.7f}"
            write_xyz(f, new_comment, coordinates, atom_types)

    # Replace original file only after successful write
    os.replace(temp_output, xyz_file)

def split_multi_xyz(multi_xyz_path, output_dir, logger=None):
    """
    Splits a multi-XYZ file into individual XYZ files.
    """
    individual_xyz_paths = []
    structures = read_xyz(multi_xyz_path, logger)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for i, (atom_count, comment, coords, atom_types) in enumerate(structures):
        output_file_path = os.path.join(output_dir, f"structure_{i+1:04d}.xyz")
        with open(output_file_path, 'w') as f:
            write_xyz(f, comment, coords, atom_types)
        individual_xyz_paths.append(output_file_path)

    if logger:
        logger.info(f"Split {multi_xyz_path} into {len(individual_xyz_paths)} individual XYZ files in {output_dir}")
    return individual_xyz_paths

def split_multi_pdb(multi_pdb_path, output_dir, logger=None):
    """
    Splits a multi-PDB file into individual PDB files, one for each model.
    Returns a list of paths to the individual PDB files.
    """
    individual_pdb_paths = []
    current_model_lines = []
    model_count = 0
    
    try:
        with open(multi_pdb_path, 'r') as f:
            for line in f:
                if line.startswith("MODEL"):
                    # If this is not the very first model, save the previous one
                    if current_model_lines:
                        model_count += 1
                        output_file_path = os.path.join(output_dir, f"model_{model_count:06d}.pdb")
                        with open(output_file_path, 'w') as out_f:
                            out_f.writelines(current_model_lines)
                        individual_pdb_paths.append(output_file_path)
                        current_model_lines = []
                    current_model_lines.append(line) # Start new model
                elif line.startswith("ENDMDL"):
                    current_model_lines.append(line)
                    model_count += 1
                    output_file_path = os.path.join(output_dir, f"model_{model_count:06d}.pdb")
                    with open(output_file_path, 'w') as out_f:
                        out_f.writelines(current_model_lines)
                    individual_pdb_paths.append(output_file_path)
                    current_model_lines = []
                else:
                    current_model_lines.append(line)
            
            # Save any remaining lines if the file doesn't end with ENDMDL after the last MODEL
            if current_model_lines and not individual_pdb_paths: # Case for single PDB without MODEL/ENDMDL
                 model_count += 1
                 output_file_path = os.path.join(output_dir, f"model_{model_count:06d}.pdb")
                 with open(output_file_path, 'w') as out_f:
                     out_f.writelines(current_model_lines)
                 individual_pdb_paths.append(output_file_path)
            elif current_model_lines and individual_pdb_paths and not individual_pdb_paths[-1].endswith(f"model_{model_count:06d}.pdb"):
                 # This handles cases where the last model doesn't have an ENDMDL
                 model_count += 1
                 output_file_path = os.path.join(output_dir, f"model_{model_count:06d}.pdb")
                 with open(output_file_path, 'w') as out_f:
                     out_f.writelines(current_model_lines)
                 individual_pdb_paths.append(output_file_path)

    except FileNotFoundError:
        if logger:
            logger.error(f"Multi-PDB file not found: {multi_pdb_path}")
        return []
    except Exception as e:
        if logger:
            logger.error(f"Error splitting multi-PDB file {multi_pdb_path}: {e}")
        return []

    if logger:
        logger.info(f"Split {multi_pdb_path} into {len(individual_pdb_paths)} individual PDB files.")
    return individual_pdb_paths

def write_multi_pdb(pdb_paths, output_path):
    """
    Merge a list of individual PDB files into a single multi-model PDB file.

    This function reads each input PDB file, extracts its "HETATM" lines until a "TER" line
    (if present), and writes them to the output file as a separate model, enclosed between
    "MODEL" and "ENDMDL" lines.

    Args:
        pdb_paths (list): List of file paths to the individual PDB files.
        output_path (str): Path to the output multi-model PDB file.
    """
    with open(output_path, 'w') as out_file:
        for i, pdb_path in enumerate(pdb_paths, start=1):
            out_file.write(f"MODEL        {i}\n")
            with open(pdb_path, 'r') as in_file:
                for line in in_file:
                    if line.startswith("HETATM"):
                        out_file.write(line)
                    elif line.startswith("TER"):
                        break
            out_file.write("ENDMDL\n")

def add_dummy_atom_to_xyz(xyz_coords, atom_types, dummy_atom_element, dummy_atom_coords):
    """
    Appends a dummy atom to existing XYZ coordinates and atom types.

    Args:
        xyz_coords (np.ndarray): NumPy array of coordinates (shape: N, 3).
        atom_types (list): List of atom type strings (length: N).
        dummy_atom_element (str): Element string for the dummy atom (e.g., "XG").
        dummy_atom_coords (np.ndarray or list): Coordinates of the dummy atom (shape: 1, 3 or 3,).

    Returns:
        tuple: (modified_xyz_coords, modified_atom_types, dummy_atom_index)
               modified_xyz_coords (np.ndarray): Updated coordinates array.
               modified_atom_types (list): Updated list of atom types.
               dummy_atom_index (int): 0-based index of the newly added dummy atom.
    """
    # Ensure dummy_atom_coords is a NumPy array with the correct shape (1, 3) for vstack
    da_coords_np = np.array(dummy_atom_coords).reshape(1, 3)

    # Append to coordinates
    updated_xyz_coords = np.vstack([xyz_coords, da_coords_np])
    
    # Append to atom types
    updated_atom_types = list(atom_types) # Ensure it's a mutable list copy
    updated_atom_types.append(dummy_atom_element)
    
    # The index of the new atom is the last index of the updated list
    dummy_atom_index = len(updated_atom_types) - 1
    
    return updated_xyz_coords, updated_atom_types, dummy_atom_index