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

def convert_optimised_arr_xyz_to_pdb(output_pdb_path, xyz_data, host_atom_count, guest_names_ordered, ingredient_map, logger=None):
    """
    Writes a PDB file from XYZ data, assigning residue names based on host and ordered guests.

    Args:
        output_pdb_path (str): Path to the output PDB file.
        xyz_data (tuple): A single (atom_count, comment, coordinates, atom_types) tuple
                          from read_xyz for the combined pull.xyz file.
        host_atom_count (int): Number of atoms in the host molecule.
        guest_names_ordered (list): Ordered list of guest names (e.g., ['ser', 'arg'])
                                    as extracted from the arrangement_N.xyz comment.
        ingredient_map (dict): Map of ingredient names (str) to Ingredient objects.
        logger (logging.Logger, optional): Logger instance. Defaults to None.
    """
    total_xyz_atoms, _, coords, atom_types = xyz_data

    # Calculate expected total atoms from host and known guests for verification
    expected_total_atoms_sum = host_atom_count
    for guest_name in guest_names_ordered:
        guest_obj = ingredient_map.get(guest_name)
        if guest_obj:
            # Call get_atom_count from imported utils module
            guest_count = get_atom_count(guest_obj.path, logger=logger)
            if guest_count == 0:
                logger.warning(f"Could not get atom count for guest '{guest_name}' from its path '{guest_obj.path}'. "
                               f"This guest's atoms will not be correctly assigned if its file is empty or invalid.")
            expected_total_atoms_sum += guest_count
        else:
            logger.warning(f"Guest '{guest_name}' from comment not found in ingredient_map. "
                           f"Its atoms will not be correctly assigned in PDB. Check config.yaml.")
            # If guest not in map, we can't get its atom count, so we just assume 0 for sum verification.

    if expected_total_atoms_sum != total_xyz_atoms:
        logger.warning(f"Total atom count mismatch for {output_pdb_path}: "
                       f"Expected sum from host/guests ({expected_total_atoms_sum}) != "
                       f"atoms in pull.xyz ({total_xyz_atoms}). "
                       f"Atom assignments in PDB might be incorrect.")

    current_atom_global_idx = 0
    current_residue_serial = 1 # PDB residue serial number
    
    with open(output_pdb_path, 'w') as f:
        # Write HOST atoms
        host_res_name = "UNK" # Standard residue name for the host
        for i in range(host_atom_count):
            if current_atom_global_idx >= total_xyz_atoms:
                logger.error(f"Ran out of atoms in XYZ data while writing host. PDB might be incomplete.")
                break
            
            atom_element = atom_types[current_atom_global_idx]
            x, y, z = coords[current_atom_global_idx]
            
            # PDB ATOM/HETATM record format (simplified example):
            # HETATM  AtomIdx AtomName ResName ChainID ResSeq   X      Y      Z      Occ  TempFactor Element
            # HETATM     1  C           HST A    1      -1.234   0.567   2.890  1.00  0.00           C
            line = (f"HETATM{current_atom_global_idx + 1:5d}  " # Atom serial number
                    f"{atom_element:<4s}" # Atom name (element symbol, left-justified)
                    f"{host_res_name:<3s} " # Residue name (3 chars, left-justified)
                    f"A" # Chain ID (fixed to A for simplicity)
                    f"{current_residue_serial:4d}    " # Residue sequence number
                    f"{x:8.3f}{y:8.3f}{z:8.3f}" # X, Y, Z coordinates (8 chars, 3 decimal places)
                    f"  1.00  0.00          " # Occupancy, Temp Factor
                    f"{atom_element:>2s}\n") # Element symbol (right-justified)
            f.write(line)
            current_atom_global_idx += 1
        current_residue_serial += 1 # Increment residue number for the next molecule

        # Write GUEST atoms
        for guest_name in guest_names_ordered:
            guest_obj = ingredient_map.get(guest_name)
            if guest_obj:
                # Call get_atom_count from imported utils module
                guest_atom_count = get_atom_count(guest_obj.path, logger=logger)
                if guest_atom_count == 0:
                    logger.warning(f"Guest '{guest_name}' has 0 atoms according to its file. Skipping in PDB.")
                    continue

                # Use uppercase first 3 characters of guest name for PDB residue name
                guest_res_name = guest_name.upper()[:3] 
                
                for i in range(guest_atom_count):
                    if current_atom_global_idx >= total_xyz_atoms:
                        logger.error(f"Ran out of atoms in XYZ data while writing guest '{guest_name}'. PDB might be incomplete.")
                        break

                    atom_element = atom_types[current_atom_global_idx]
                    x, y, z = coords[current_atom_global_idx]
                    
                    if i == 1 and (atom_element.upper() == "C" or atom_element.upper() == "CA" or atom_element.upper() == "CX"):
                        atom_element = "CA"
                    line = (f"ATOM  {current_atom_global_idx + 1:5d}  "
                            f"{atom_element:<4s}"
                            f"{guest_res_name:<3s} "
                            f"A"
                            f"{current_residue_serial:4d}    "
                            f"{x:8.3f}{y:8.3f}{z:8.3f}"
                            f"  1.00  0.00          "
                            f"{atom_element[0]:>2s}\n")
                    f.write(line)
                    current_atom_global_idx += 1
                current_residue_serial += 1 # Increment residue number for the next guest
            else:
                logger.error(f"Ingredient '{guest_name}' not found in ingredient_map. Cannot write its atoms to PDB.")
                # If guest not found, we cannot determine its atom count or path, so skip its PDB writing.

        f.write("END\n") # Standard PDB file terminator
        logger.info(f"Successfully wrote PDB file: {output_pdb_path}")
        
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
    Splits a multi-structure PDB file into individual PDB files,
    including the CONECT records from the end of the multi-structure file
    in each split file.

    Args:
        multi_pdb_path (str): Path to the input multi-structure PDB file.
        output_dir (str): Directory to save the split PDB files.
        logger (logging.Logger, optional): Logger instance. Defaults to None.

    Returns:
        list: A list of paths to the individual PDB files created.
    """
    individual_pdb_paths = []
    current_model_lines = []
    model_count = 0
    conect_section = []

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        if logger:
            logger.info(f"Created output directory: {output_dir}")

    try:
        with open(multi_pdb_path, 'r') as f:
            lines = f.readlines()

        # --- Extract CONECT records from the end of the file ---
        # Find the start of CONECT records by searching backwards
        conect_records_start_index = -1
        # Iterate backwards to find the block of CONECT records
        for i in range(len(lines) - 1, -1, -1):
            if lines[i].startswith("CONECT"):
                conect_records_start_index = i
            elif conect_records_start_index != -1 and not lines[i].strip():
                # Found an empty line after a CONECT, so the CONECT block starts from conect_records_start_index
                break
            elif conect_records_start_index != -1 and not lines[i].startswith("CONECT"):
                 # Found a non-CONECT line after a CONECT, means the CONECT block started earlier
                 # This handles cases where there might be non-CONECT lines between CONECT blocks or before the first CONECT.
                 # We want the *last contiguous block* of CONECT records.
                 conect_records_start_index = i + 1 # The block starts from the next line
                 break


        if conect_records_start_index != -1:
            conect_section = [line for line in lines[conect_records_start_index:] if line.startswith("CONECT")]
            # Remove CONECT records from the main lines list to avoid re-processing them
            # We only remove them if they were found at the very end of the file as expected.
            if conect_records_start_index == len(lines) - len(conect_section): # Check if it's the last block
                lines = lines[:conect_records_start_index]
        else:
            if logger:
                logger.warning("No CONECT records found at the end of the PDB file.")

        # --- Process models ---
        for line in lines:
            if line.startswith("MODEL"):
                # If we were already processing a model, save it before starting a new one
                if current_model_lines:
                    model_count += 1
                    output_file_name = f"model_{model_count:06d}.pdb"
                    output_file_path = os.path.join(output_dir, output_file_name)
                    
                    with open(output_file_path, 'w') as out_file:
                        out_file.writelines(current_model_lines)
                        # Add CONECT records after the model data and its ENDMDL
                        for conect_line in conect_section:
                            out_file.write(conect_line)
                    individual_pdb_paths.append(output_file_path)
                    current_model_lines = [] # Reset for the next model
                current_model_lines.append(line) # Add the current MODEL line to the new set
            elif line.startswith("ENDMDL"):
                current_model_lines.append(line)
                # This ENDMDL signals the end of a model, so we save it immediately
                model_count += 1
                output_file_name = f"model_{model_count:06d}.pdb"
                output_file_path = os.path.join(output_dir, output_file_name)

                with open(output_file_path, 'w') as out_file:
                    out_file.writelines(current_model_lines)
                    # Add CONECT records after the model data and its ENDMDL
                    for conect_line in conect_section:
                        out_file.write(conect_line)
                individual_pdb_paths.append(output_file_path)
                current_model_lines = [] # Reset for the next model
            else:
                # Accumulate lines that are part of the current model (ATOM, HETATM, REMARK, etc.)
                current_model_lines.append(line)

        # Handle the case where the file might not end with ENDMDL after the last MODEL
        # or if it's a single model file without MODEL/ENDMDL
        if current_model_lines:
            # Check if this model has already been saved (i.e., if it ended with ENDMDL)
            # This logic can be tricky if the file format is inconsistent.
            # A more robust check: if the last saved path corresponds to the current model_count
            # then it's already saved.
            
            # Simple check for now: if current_model_lines isn't empty, and it wasn't
            # explicitly ended by ENDMDL (which would have cleared it), save it.
            # This primarily catches the very last model if it lacks an ENDMDL or if it's a single model file.
            if not individual_pdb_paths or (model_count == 0 and current_model_lines): # Case: single PDB or last model of multi-PDB
                model_count += 1
                output_file_name = f"model_{model_count:06d}.pdb"
                output_file_path = os.path.join(output_dir, output_file_name)
                
                with open(output_file_path, 'w') as out_file:
                    out_file.writelines(current_model_lines)
                    # Add CONECT records
                    for conect_line in conect_section:
                        out_file.write(conect_line)
                    # Ensure ENDMDL if not present for single-model files or last model
                    if not any(line.startswith("ENDMDL") for line in current_model_lines):
                        out_file.write("ENDMDL\n")
                individual_pdb_paths.append(output_file_path)

    except FileNotFoundError:
        logger.error(f"Multi-PDB file not found: {multi_pdb_path}")
        return []
    except Exception as e:
        logger.error(f"Error splitting multi-PDB file {multi_pdb_path}: {e}")
        return []

    logger.info(f"Split {multi_pdb_path} into {len(individual_pdb_paths)} individual PDB files, each with CONECT records.")
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
                    if line.startswith("HETATM") or line.startswith("ATOM"):
                        out_file.write(line.replace("ATOM  ", "HETATM"))
                    elif line.startswith("TER"):
                        break
            out_file.write("ENDMDL\n")

def add_dummy_atom_to_xyz(xyz_coords, atom_types, dummy_atom_coords):
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
    updated_atom_types.append("DA")
    
    # The index of the new atom is the last index of the updated list
    dummy_atom_index = len(updated_atom_types) - 1
    
    return updated_xyz_coords, updated_atom_types, dummy_atom_index