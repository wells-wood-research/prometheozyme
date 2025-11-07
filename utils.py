import numpy as np
import io
import os
import pandas as pd
import re
import pdbUtils
from typing import Optional, Tuple, Any
import logging

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
    if not pathXYZ.lower().endswith(".pdb"):
        print(f"File must be XYZ type (.xyz), provided file is: {pathXYZ}")
        return False
    if not os.path.exists(pathXYZ):
        print(f"File does not exist: {pathXYZ}")
        return False
    return True

def pdb2df(pathPDB: str, logger: Optional[logging.Logger] = None) -> Optional[pd.DataFrame]:
    """Convert a PDB file to a DataFrame."""
    loggger = get_logger()

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


def df2xyz(df: pd.DataFrame, comment: str, outXYZ: str) -> None:
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
        df2xyz(df, comment=structure_name, outfile=outXYZ)
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
                structures.append((atom_count, comment, coordinates, atom_types)) # TODO is that really the best name? Why is it a list?
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

def split_docker_results(multi_xyz_path, logger=None):
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

    for i, (atom_count, comment, coords, atom_types) in enumerate(structures, start=1):
        # TODO use atom_count to check if it agrees with the total of ingredients atom count?
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