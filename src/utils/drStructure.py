import io
import os
import re
import logging
import numpy as np
import pandas as pd
from pdbUtils import pdbUtils
from typing import Tuple, Optional, Any

from utils.drEval import evaluate_distance

########################
## LOGGING
########################

def get_logger(logger: Optional[logging.Logger] = None) -> logging.Logger:
    """Return a valid logger, creating one if not provided."""
    if logger is None:
        logger = logging.getLogger(__name__)
        if not logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')
            handler.setFormatter(formatter)
            logger.addHandler(handler)
        logger.setLevel(logging.INFO)
    return logger

########################
## FILE TYPE CONVERSION
########################

def isPDB(pathPDB: str) -> bool:
    if not isinstance(pathPDB, str):
        print("pathPDB must be a string.")
        return False
    if not pathPDB.lower().endswith(".pdb"):
        print(f"File must be PDB type (.pdb), provided file is: {pathPDB}")
        return False
    if not os.path.exists(pathPDB):
        print(f"File does not exist: {pathPDB}")
        return False
    return True

def isXYZ(pathXYZ: str) -> bool:
    if not isinstance(pathXYZ, str):
        print("pathXYZ must be a string.")
        return False
    if not pathXYZ.lower().endswith(".xyz"):
        print(f"File must be XYZ type (.xyz), provided file is: {pathXYZ}")
        return False
    if not os.path.exists(pathXYZ):
        print(f"File does not exist: {pathXYZ}")
        return False
    return True

def pdb2df(pathPDB: str, logger: Optional[logging.Logger] = None) -> Optional[pd.DataFrame]:
    """Convert a PDB file to a DataFrame."""
    logger = get_logger()

    if not isPDB(pathPDB):
        return None

    try:
        df = pdbUtils.pdb2df(pathPDB)
        if not isinstance(df, pd.DataFrame):
            logger.error("pdbUtils.pdb2df did not return a DataFrame.")
            return None
        return df
    except Exception as e:
        logger.exception(f"Exception during PDB → DataFrame conversion: {e}")
        return None


def df2pdb(df: pd.DataFrame, outPDB: str, logger: Optional[logging.Logger] = None) -> Optional[Any]:
    """Convert a DataFrame to a PDB file."""
    logger = get_logger(logger)

    if not isinstance(df, pd.DataFrame) or df.empty:
        logger.error("Invalid or empty DataFrame provided.")
        return None

    try:
        pdbUtils.df2pdb(df=df, outFile=outPDB)
        if not isPDB(outPDB):
            logger.error("pdbUtils.df2pdb did not return a correct PDB file.")
            return None
    except Exception as e:
        logger.exception(f"Exception during DataFrame → PDB conversion: {e}")
        return None


def df2xyz(df: pd.DataFrame, comment: str, outXYZ: str, logger: Optional[logging.Logger] = None) -> None:
    """Write DataFrame to XYZ file format."""
    logger = get_logger(logger)

    if not isinstance(df, pd.DataFrame) or df.empty:
        raise ValueError("Invalid or empty DataFrame provided.")

    if not all(col in df.columns for col in ["ELEMENT", "X", "Y", "Z"]):
        raise ValueError("DataFrame must contain columns: ELEMENT, X, Y, Z")

    with open(outXYZ, "w") as f:
        f.write(f"{len(df)}\n")
        f.write(f"{comment or ''}\n")
        for _, row in df.iterrows():
            f.write(f"{row['ELEMENT']} {row['X']:27.17f} {row['Y']:27.17f} {row['Z']:27.17f}\n")
    
    if not isXYZ(outXYZ):
        logger.error("df2xyz did not return a correct XYZ file.")
        return None

def pdb2xyz(pathPDB: str, outXYZ: str, logger: Optional[logging.Logger] = None) -> None:
    """Convert a PDB file to XYZ format."""
    logger = get_logger(logger)

    df = pdb2df(pathPDB, logger)
    if df is None:
        logger.error("Failed to load PDB file; aborting conversion.")
        return

    try:
        structure_name = os.path.basename(pathPDB).split(".")[0]
        df2xyz(df, comment=structure_name, outXYZ=outXYZ)
        if not isXYZ(outXYZ):
            logger.error("pdb2xyz did not return a correct XYZ file.")
            return None
    except Exception as e:
        logger.exception(f"Exception during PDB → XYZ conversion: {e}")


def xyz2df(pathXYZ: str, logger: Optional[logging.Logger] = None) -> Optional[Tuple[int, str, pd.DataFrame]]:
    """Convert an XYZ file to a DataFrame."""
    logger = get_logger(logger)

    if not isXYZ(pathXYZ):
        return None

    try:
        atom_count, comment, coords, atom_types = read_xyz(pathXYZ, logger)
        df = pd.DataFrame({
            "ELEMENT": atom_types,
            "X": coords[:, 0],
            "Y": coords[:, 1],
            "Z": coords[:, 2],
        })
        return atom_count, comment, df
    except Exception as e:
        logger.exception(f"Exception during XYZ → DataFrame conversion: {e}")
        return None


def xyz2pdb(pathXYZ: str, outPDB: str, logger: Optional[logging.Logger] = None) -> None:
    """Convert an XYZ file to a PDB format."""
    logger = get_logger(logger)

    xyz_data = xyz2df(pathXYZ, logger)
    if xyz_data is None:
        logger.error("Failed to read XYZ file; aborting conversion.")
        return

    _, _, df = xyz_data
    try:
        df2pdb(df, outPDB, logger)
        if not isXYZ(outPDB):
            logger.error("xyz2pdb did not return a correct PDB file.")
            return None
    except Exception as e:
        logger.exception(f"Exception during XYZ → PDB conversion: {e}")

########################
## FILE READING
########################

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

def parse_energy_comment(comment):
    eopt = einter = None
    if comment:
        eopt_match = re.search(r"Eopt=(-?\d+\.\d+)", comment)
        einter_match = re.search(r"Einter=(-?\d+\.\d+)", comment)
        if eopt_match: eopt = float(eopt_match.group(1))
        if einter_match: einter = float(einter_match.group(1))
    return eopt, einter

def extract_ok_docker_results(multi_xyz_path, n_atoms_host, biases, logger=None):
    """
    - Split a multi-structure XYZ file.
    - Extract Eopt/Einter from the comment line.
    - Save each as a single XYZ file.
    - Return list of (Result, DataFrame).
    """
    base_dir = os.path.dirname(multi_xyz_path)
    base_name = os.path.basename(multi_xyz_path)

    # Reuse universal parser
    structures = read_xyz(multi_xyz_path, logger=logger)

    results = []

    for i, (atom_count, comment, coords, atom_types) in enumerate(structures):
        # TODO after optmisation there might be duplicated results - need to remove based on RMSD
        # use AMPAL? https://isambard-uob.github.io/ampal/ampal.html#ampal.base_ampal.BaseAmpal.rmsd
        allOk = True
        for bias in biases:
            [guestIdx, hostIdx] = bias.get("atoms", [0, 0])
            val = bias.get("val", 0)
            tol = bias.get("tol", 0.5)
            distance = evaluate_distance(coords, guestIdx, hostIdx, n_atoms_host)
            isOk = val - tol <= distance <= val + tol
            logger.info(f"Distance between restrained atoms: {distance}, expected {val} within {tol} tolerance ({'failed' if not isOk else 'passed'}).")
            if not isOk:
                allOk = False
                break
        if not allOk:
            continue

        # --- Parse energies using regex ---
        eopt = einter = None
        if comment:
            eopt, einter = parse_energy_comment(comment)

        # --- Create dataframe ---
        df = pd.DataFrame({
            "ELEMENT": atom_types,
            "X": coords[:, 0],
            "Y": coords[:, 1],
            "Z": coords[:, 2],
        })

        # --- Construct new output filename ---
        new_name = re.sub(r"\.all\.optimized\.xyz$", "", base_name)
        new_path = os.path.join(base_dir, f"{new_name}.result.{i}.xyz")

        # --- Reuse write_xyz() ---
        write_xyz(new_path, comment, coords, atom_types)

        results.append((new_path, eopt, einter, df))

    if logger:
        logger.info(f"Split {multi_xyz_path} into {len(results)} results.")

    if len(results) < 1:
        logger.warning("No docker result structures passed the restraint test for this host seed. This is not critical unless no subsequent results exist. To increase success rate, consider increasing tolerance or force on the relevant restraint.")
    
    return results

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
                    if line.startswith(("ATOM", "HETATM")):
                        out_file.write(line)
                    out_file.write(line)
            out_file.write("\nENDMDL\n")
        out_file.write("END\n")

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