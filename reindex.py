import pdbUtils
import pandas as pd
import argparse
import os
import re
from typing import Optional, List

from rdkit import Chem as RDchem
from rdkit.Chem import Draw
from rdkit.Chem import rdDetermineBonds
from rdkit.Chem import AllChem as RDchem
from rdkit.Chem import rdFMCS
from rdkit.Chem import rdForceFieldHelpers as RDforceFields
from rdkit.Chem import Descriptors

class XYZFileFormatError(Exception):
    """Exception raised for errors in the XYZ file format."""
    pass

class SubstructureNotFound(Exception):
    """Exception raised when RDKit's findMCS (maximum common substructure) calculation failed to finish."""
    pass

class StructureNotOptimised(Exception):
    """Structure optimisation with rdkit.Chem.rdForceFieldHelpers.UFFOptimizeMolecule did not converge"""

def load_molecule(molec: str) -> tuple[pd.DataFrame, RDchem.Mol]:
    """
    Read molecule structure file into a pandas dataframe, handing PDB and XYZ file formats.

    Args:
        molec (str): Absolute path to a molecule structure file.

    Returns:
        pd.DataFrame:   Dataframe describing the input structure with columns:
                        ['ATOM', 'ATOM_ID', 'ATOM_NAME', 'RES_NAME', 'CHAIN_ID', 'RES_ID', 'X', 'Y', 'Z', 'OCCUPANCY', 'BETAFACTOR', 'ELEMENT']
        RDchem.Mol:     RDKit molecule object constructred from the input structure file
    """
    name, extension = get_file_format(molec)
    if extension == "pdb":
        df = pdbUtils.pdb2df(molec)
        mol = RDchem.rdmolfiles.MolFromPDBFile(molec)
    else:
        # Assume .xyz as checked at get_file_format level
        df = xyz2df(molec)
        mol = RDchem.rdmolfiles.MolFromXYZFile(molec)
        RDchem.rdDetermineBonds.DetermineBonds(mol,charge=1) # TODO Must allow user to input their charge
    mol.SetProp("name", name)
    return df, mol, name

def get_file_format(file: str) -> Optional[tuple[str, str]]:
    """
    Check file format by matching extension string against implemented options.

    Args:
        file (str): Absolute path to a molecule structure file.

    Returns:
        Optional[str]: Extension format of the input file without the leading dot,
                       or None if the format is not supported.

    Raises:
        TypeError: If the input is not a string.
        ValueError: If the file path is empty or has no extension.
    """
    # Type checking
    if not isinstance(file, str):
        raise TypeError("Input must be a string representing a file path.")

    # Check if file path is empty
    if not file.strip():
        raise ValueError("File path cannot be empty.")

    if not is_valid_path(file):
        raise ValueError("File path cannot be empty.")
    try:
        # Get the extension (with the dot)
        root, extension = os.path.splitext(file)
        name = os.path.basename(root)
        # Check if there is no extension
        if not extension:
            raise ValueError("Please provide a full file path with an extension (e.g., filename.pdb or filename.xyz).")

        # Define supported extensions (as lowercase for consistency)
        global IMPLEMENTED_EXTENSIONS
        IMPLEMENTED_EXTENSIONS = {".pdb", ".xyz"}  # Using a set for faster lookup
        if extension not in IMPLEMENTED_EXTENSIONS:
            raise ValueError(f"Unsupported file extension '{extension}'. Allowed extensions: {sorted(IMPLEMENTED_EXTENSIONS)}")

        return name, extension.lstrip('.')

    except AttributeError:
        raise TypeError("File path must be a string or path-like object.")
    except Exception as e:
        raise RuntimeError(f"Error processing file: {str(e)}") from e

def is_valid_path(file: str) -> bool:
    '''
    Check if file exists at the provided path.

    Arguments:
        file: absolute path to molecule structure file
    Returns:
        bool: True if file exists, False otherwise
    '''
    exists = os.path.isfile(file)
    return exists

def xyz2df(molec: str) -> pd.DataFrame:
    """
    CREDIT: Inspired by pdbUtils (https://github.com/ESPhoenix/pdbUtils) code by Eugene Shrimpton-Pheonix.

    Read .xyz file into a pandas dataframe.

    Args:
        molec (str): Absolute path to a molecule structure file.

    Returns:
        pd.DataFrame:   Dataframe describing the input structure with columns:
                        ['ATOM', 'ATOM_ID', 'ATOM_NAME', 'RES_NAME', 'CHAIN_ID', 'RES_ID', 'X', 'Y', 'Z', 'OCCUPANCY', 'BETAFACTOR', 'ELEMENT']

    Raises:
        FileNotFoundError: If the XYZ file does not exist.
        XYZFileFormatError: If the XYZ file format is invalid (e.g., missing header lines).
        RuntimeError: If there are issues reading the file or parsing data.

    Example:
        >>> xyz2df("/home/mchrnwsk/reindexer/src/tests/water.xyz")
            ATOM  ATOM_ID ATOM_NAME RES_NAME CHAIN_ID  RES_ID      X       Y    Z  OCCUPANCY  BETAFACTOR ELEMENT
        0  HETATM        1         O      UNK        A       1  0.000 -0.0589  0.0        1.0         0.0       O
        1  HETATM        2         H      UNK        A       2 -0.811  0.4677  0.0        1.0         0.0       H
        2  HETATM        3         H      UNK        A       3  0.811  0.4677  0.0        1.0         0.0       H
    """
    columns = ['ATOM', 'ATOM_ID', 'ATOM_NAME', 'RES_NAME', 'CHAIN_ID', 'RES_ID', 'X', 'Y', 'Z', 'OCCUPANCY', 'BETAFACTOR', 'ELEMENT']
    data: List[List] = []

    try:
        with open(molec, 'r') as xyz_file:
            lines = xyz_file.readlines()

            # Check that first two lines are for number of atoms and comment
            if len(lines) < 2:
                raise XYZFileFormatError("Please provide .xyz files with header lines: total number of atoms in line 1, and description in line 2.")

            # First line should be the number of atoms
            try:
                noAtms = int(lines[0].strip())
            except ValueError:
                raise XYZFileFormatError("First line of XYZ file must be an integer (number of atoms).")

            # Verify that the number of remaining lines matches the number of atoms
            atom_lines = lines[2:]  # Skip the first two header lines
            if len(atom_lines) != noAtms:
                raise XYZFileFormatError(f"Expected {noAtms} atoms, but found {len(atom_lines)} data lines.")

            # Process atom data starting from line 3
            for line_num, line in enumerate(atom_lines, start=3):  # Start numbering from 3
                line = line.strip()
                if not line:  # Skip empty lines
                    continue

                # Split the line into columns (assuming space-separated)
                parts = line.split()
                if len(parts) < 4:  # Ensure at least atom name, x, y, z are present
                    raise XYZFileFormatError(f"Line {line_num} in XYZ file has insufficient columns. Expected at least 4 (atom, x, y, z).")

                atom_name = parts[0].strip()  # First column: atom name or element
                try:
                    x = float(parts[1].strip())  # Second column: x coordinate
                    y = float(parts[2].strip())  # Third column: y coordinate
                    z = float(parts[3].strip())  # Fourth column: z coordinate
                except (ValueError, IndexError) as e:
                    raise XYZFileFormatError(f"Invalid coordinate data on line {line_num}: {str(e)}")

                # Fill other columns as specified
                atom_type = "HETATM"
                res_name = "UNK"
                chain_id = "A"
                atom_id = line_num - 2  # Line number minus header lines (1-based for users)
                res_id = line_num - 2   # Same as atom_id for simplicity
                occupancy = 1.00
                temp_factor = 0.00
                element = atom_name  # Element is same as atom name

                data.append([atom_type, atom_id, atom_name, res_name, chain_id, res_id, x, y, z, occupancy, temp_factor, element])

        if not data:
            raise XYZFileFormatError("No valid atom data found in XYZ file.")

        return pd.DataFrame(data, columns=columns)

    except FileNotFoundError:
        raise FileNotFoundError(f"XYZ file not found: {molec}")
    except Exception as e:
        raise RuntimeError(f"Error reading XYZ file: {str(e)}") from e

def df2xyz(df: pd.DataFrame, filepath: str) -> None:
    """
    Save a pandas DataFrame containing atomic structure data into an XYZ file.

    The first line contains the number of atoms (rows in the DataFrame), the second line is a comment line, and subsequent lines 
    list the atom label (from the 'ELEMENT' column of the dataframe) and X, Y, Z coordinates, separated by tabs.

    Args:
        df (pd.DataFrame):  DataFrame describing the atomic structure with columns:
                            ['ATOM', 'ATOM_ID', 'ATOM_NAME', 'RES_NAME', 'CHAIN_ID', 'RES_ID', 'X', 'Y', 'Z', 
                            'OCCUPANCY', 'BETAFACTOR', 'ELEMENT']. Must include at least 'ELEMENT', 'X', 'Y', 'Z'.
        filepath (str):     Absolute path where the XYZ file will be saved, including the '.xyz' extension.

    Raises:
        KeyError: If required columns ('ELEMENT', 'X', 'Y', 'Z') are missing from the DataFrame.
        IOError: If there is an issue writing to the specified file path.
    """
    # Check for required columns
    colRequired = ['ELEMENT', 'X', 'Y', 'Z']
    if not all(col in df.columns for col in colRequired):
        raise KeyError(f"DataFrame is missing required columns: {colRequired}")

    # Number of atoms is the number of rows in the DataFrame
    noAtms = len(df)

    # Prepare the header lines
    header = f"{noAtms}\nSaved from a dataframe by https://github.com/wells-wood-research/reindexer.git\n"

    # Format each row as "ELEMENT    X    Y    Z" with tab separation and fixed-width floating-point notation
    atom_lines = []
    for _, row in df.iterrows():
        element = f"{row['ELEMENT']:<2}"  # Left-align element, fixed width of 2 (e.g., 'C ', 'N ')
        x = f"{row['X']:>9.5f}"          # Right-align X, width of 9, 5 decimal places
        y = f"{row['Y']:>9.5f}"          # Right-align Y, width of 9, 5 decimal places
        z = f"{row['Z']:>9.5f}"          # Right-align Z, width of 9, 5 decimal places
        line = f"{element}\t{x}\t{y}\t{z}"
        atom_lines.append(line)

    # Combine header and atom lines
    xyz_content = header + "\n".join(atom_lines) + "\n"  # Add final newline for consistency

    # Write to file
    try:
        with open(filepath, 'w') as xyz_file:
            xyz_file.write(xyz_content)
    except IOError as e:
        raise IOError(f"Error writing to XYZ file at {filepath}: {str(e)}")

#############################################################################################################

def get_maximum_common_substructure(mol1: RDchem.Mol, mol2: RDchem.Mol) -> RDchem.Mol:
    """
    CREDIT: from RDKit cookbook (https://www.rdkit.org/docs/Cookbook.html#highlight-molecule-differences) code by Takayuki Serizawa.
    
    Finds maximum common substructure between two molecules:
    the largest substructure common to both input molecules.
    FindMCS returns MCSResult object, including the following properties:
        canceled:   if True, the MCS calculation did not finish;    ==> abort
        queryMol:   query molecule for the MCS;                     ==> return object

    Args:
        mol1 (RDchem.Mol):      The reference RDKit molecule object.
        mol2 (RDchem.Mol):      The referee RDKit molecule object to compare against mol1.

    Returns:
        RDchem.Mol:             The maximum common substructure as an RDKit molecule object.

    Raises:
        SubstructureNotFound: If MCS cannot be found by RDKit.
    """

    mcs_params = rdFMCS.MCSParameters()
    mcs_params.AtomTyper = rdFMCS.AtomCompare.CompareAny  # Match atoms by connectivity only
    mcs_params.BondTyper = rdFMCS.BondCompare.CompareAny  # Match bonds by connectivity only
    mcs_params.BondCompareParameters.MatchStereo = False  # Explicitly ignore stereochemistry
    mcs_params.AtomCompareParameters.MatchChiralTag = False  # Ignore chirality (already default)
    mcs_params.Timeout = 3600  # Set timeout (in seconds)
    mcs_params.Verbose = False  # Enable verbose output for debugging

    mcs = rdFMCS.FindMCS([mol1,mol2], mcs_params)
    if mcs.canceled:
        raise SubstructureNotFound("Common substructure could not be found. Investigate your inputs.")
    mcs_mol = mcs.queryMol
    return mcs_mol

#############################################################################################################

def save_image_difference(mol1: RDchem.Mol, match1: tuple[int], mol2: RDchem.Mol, match2: tuple[int], outpath: str) -> None:
    """
    Save an SVG image comparing two molecules, visually highlighting their differences.

    This function generates and saves an SVG image showing the reference and referee 
    molecules side by side, with atom indices labeled and differences highlighted in red. 
    The image is resized to fit the content without excessive margins.

    Args:
        mol1 (RDchem.Mol):      The reference molecule object (expected to be RDchem.Mol).
        match1 (tuple[int]):    The atom indices in mol1 that match the MCS.
        mol2 (RDchem.Mol):      The referee molecule object (expected to be RDchem.Mol).
        match2 (tuple[int]):    The atom indices in mol2 that match the MCS.
        mcs_mol (RDchem.Mol):   The maximum common substructure molecule (RDchem.Mol).
        outpath (str):          The output path (without extension) where the SVG file will be saved.

    Returns:
        None
    """
    label_atom_indices(mol1)
    label_atom_indices(mol2)
    diff1, diff2 = get_atom_difference(mol1, match1), get_atom_difference(mol2, match2)

    mols = [mol1, mol2]
    sub_width = 500
    svg_width = len(mols)*sub_width
    sub_height = 500
    svg_height = sub_height

    img = RDchem.Draw.MolsToGridImage(mols = mols
                                      , subImgSize = (sub_width, sub_height)
                                      , legends = [f"reference: {mol1.GetProp('name')}", f"referee: {mol2.GetProp('name')}"]
                                      , highlightAtomLists=[diff1, diff2]
                                      , useSVG = True
                                      )
    resized_img = edit_image_dimensions(img, svg_width, svg_height)

    with open(f"{outpath}.svg", 'w') as f:
        f.write(resized_img)

def edit_image_dimensions(img: str, width: float, height: float) -> str:
    """
    Resize an SVG image to fit its content by adjusting width, height, and viewBox.

    This function modifies the SVG string to ensure the image dimensions match the 
    content, removing unnecessary margins. It updates the 'width', 'height', and 
    'viewBox' attributes, as well as the background rectangle.

    Args:
        img (str):      input SVG image as a string.
        width (float):  width of the subimages
        height (float): height of the subimages

    Returns:
        str:            The modified SVG string with updated dimensions and viewBox to fit the content.

    Note:
        Small padding (20 pixels) is added to the width and height to ensure no content is cut off.
    """
    # Update SVG string
    resized_img = re.sub(r'width=\'[^\']+\' height=\'[^\']+\' viewBox=\'[^\']+\'', 
                    f'width=\'{width+20}px\' height=\'{height+20}px\' viewBox=\'0.0 0.0 {width} {height}\'', 
                    img, 1)

    # Update the svg rect
    resized_img = re.sub(r'<rect style=[^\>]+ width=\'[^\']+\' height=\'[^\']+\' x=\'[^\']+\' y=\'[^\']+\'>',
                    f'<rect style="opacity:1.0;fill:#FFFFFF;stroke:none" width="{width}" height="{height}" x="0" y="0">',
                    resized_img, 1)
    return resized_img

def label_atom_indices(mol: RDchem.Mol) -> None:
    """
    Label each atom in a molecule with its index.

    This function iterates over all atoms in the input molecule and assigns each 
    atom its index as a map number, which can be used for visualization or comparison.

    Args:
        mol (RDchem.Mol):    The molecule object to be labelled.

    Returns:
        None
    """
    for idx, atom in enumerate(mol.GetAtoms()):
        atom.SetAtomMapNum(idx)

def get_atom_difference(mol: RDchem.Mol, match: tuple[int]) -> list:
    """
    Identify atoms in two molecules, which differ from their maximum common substructure.

    This function compares the atoms in two molecules against their matches in the 
    maximum common substructure (MCS) to find atoms that do not match. These atoms 
    are returned as lists of indices, which can be used to highlight differences 
    (e.g., in red) in visualizations.

    Args:
        mol (RDchem.Mol):      molecule object
        match (tuple[int]):    list of atom indices in mol that match the MCS
    Returns:
       (list):  List of atom indices in a mol that do not match MCS of mol and another mol.
    """
    diff = []
    for atom in mol.GetAtoms():
        if atom.GetIdx() not in match:
            diff.append(atom.GetIdx())
    return diff

#############################################################################################################

def split_df_by_H(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Split a DataFrame describing molecular structure into two DataFrames: one for non-hydrogen atoms and one for hydrogen atoms.

    The split is determined by checking the 'ELEMENT' column, hydrogen atoms ahave 'H' in 'ELEMENT'.

    Args:
        df (pd.DataFrame):  Dataframe describing the input structure with columns:
                            ['ATOM', 'ATOM_ID', 'ATOM_NAME', 'RES_NAME', 'CHAIN_ID', 'RES_ID', 'X', 'Y', 'Z', 'OCCUPANCY', 'BETAFACTOR', 'ELEMENT']

    Returns:
        tuple[pd.DataFrame, pd.DataFrame]: A tuple containing two DataFrames:
            - First DataFrame: Rows where the atom is not hydrogen ('ELEMENT' is not 'H').
            - Second DataFrame: Rows where the atom is hydrogen ('ELEMENT' is 'H').

    Warning:
        If no rows contain hydrogen atoms (i.e., no rows where 'ATOM_NAME' or 'ELEMENT' is 'H'), a warning is printed 
        indicating that no hydrogen atoms were found in the DataFrame. The second DataFrame in this case will be empty.
    """
    # Check if either ATOM_NAME or ELEMENT is 'H' to identify hydrogen atoms
    isH = df['ATOM_NAME'].str.upper().str.match(r'^H\d*$')

    # Split DataFrame
    dfNotH = df[~isH].copy()  # Rows where atom is not hydrogen
    dfH = df[isH].copy()       # Rows where atom is hydrogen

    return dfNotH, dfH

def get_reindexed_dataframes(df1: pd.DataFrame, df2: pd.DataFrame, match1: tuple[int], match2: tuple[int], diff1: list, diff2: list) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Reindex two DataFrames to match indices between elements that are in common, and appending the unique elements at the end.
    Indices of the reference (df1) are maintained, while referee (df2) is reindexed. Hydrogen atoms are preserved by splitting
    and merging them back after reindexing non-hydrogen atoms.

    Args:
        df1 (pd.DataFrame):     The reference DataFrame whose indices should remain unchanged
        df2 (pd.DataFrame):     The referee DataFrame to be reindexed to match the reference
        match1 (tuple[int]):    Tuple of indices in df1 that correspond to common elements.
        match2 (tuple[int]):    Tuple of indices in df2 that correspond to the same common elements as match1.
        diff1 (list[int]):      List of indices of atoms unique to mol1, or None if no unique elements.
        diff2 (list[int]):      List of indices of atoms unique to mol2, or None if no unique elements.

    Returns:
        tuple[pd.DataFrame, pd.DataFrame]: A tuple containing:
            - First DataFrame: The reference DataFrame with unchanged indices for common non-hydrogen elements, unique non-hydrogen
                               elements appended, and hydrogen atoms merged back.
            - Second DataFrame: The reindexed referee DataFrame where common non-hydrogen elements match the reference indices,
                                unique non-hydrogen elements are appended, and hydrogen atoms are merged back.

    Raises:
        TypeError: If diff1 or diff2 is not a list or None.
    """
    # Validate diff1 and diff2
    if not (isinstance(diff1, list) or diff1 is None):
        raise TypeError(f"diff1 must be a list or None, got {type(diff1)}")
    if not (isinstance(diff2, list) or diff2 is None):
        raise TypeError(f"diff2 must be a list or None, got {type(diff2)}")

    # Split DataFrames into non-hydrogen and hydrogen atoms
    df1NotH, df1H = split_df_by_H(df1)
    df2NotH, df2H = split_df_by_H(df2)

#    print("df1NotH")
#    print(df1NotH)
#    print()
#    print("df1H")
#    print(df1H)
#    print()
#    print("df2NotH")
#    print(df2NotH)
#    print()
#    print("df2H")
#    print(df2H)
#    print()
    
    # Convert tuples to lists for easier manipulation
    idx1common = list(match1)
    idx2common = list(match2)

#    print("idx1common")
#    print(idx1common)
#    print()
#    print("idx2common")
#    print(idx2common)
#    print()

    # Treat None as empty list
    diff1 = diff1 or []
    diff2 = diff2 or []

    # --- Reindexing for common non-hydrogen atoms ---
    # Create a new DataFrame for df1's common atoms, maintaining their original indices
    df1common = df1NotH.iloc[idx1common].copy() if idx1common else pd.DataFrame(columns=df1NotH.columns)
#    print("df1common")
#    print(df1common)
#    print()

    # Create a temporary DataFrame for df2's common atoms, indexed by their corresponding df1 indices
    # This is the core reindexing step for common atoms.
    df2common_reindexed_data = []
    for ref_idx, referee_idx in zip(idx1common, idx2common): # Iterate through the paired common indices
        # Get the row from df2NotH using its original index label (referee_idx)
        # and set its new index to ref_idx
        row = df2NotH.iloc[referee_idx].copy()
        df2common_reindexed_data.append(row.rename(ref_idx)) # Rename the series index to the new index

#    print("df2common_reindexed_data")
#    print(df2common_reindexed_data)
#    print()

    if df2common_reindexed_data:
        df2common = pd.DataFrame(df2common_reindexed_data, columns=df2NotH.columns)
        # Sort by the new index to ensure consistent order
        df2common = df2common.sort_index()
    else:
        df2common = pd.DataFrame(columns=df2NotH.columns)

#    print("df2common")
#    print(df2common)
#    print()

    # --- Handle additional non-hydrogen elements for df1 ---
    df1_parts = [df1common]
    if diff1:
        df1diff = df1NotH.loc[diff1].copy()
        # Find the maximum index already present in df1common, or -1 if empty
        idxMaxCurrent = max(df1common.index) if not df1common.empty else -1
        # Assign new sequential indices starting from after the current max index
        new_diff1_indices = range(idxMaxCurrent + 1, idxMaxCurrent + 1 + len(diff1))
        df1diff.index = new_diff1_indices
        df1_parts.append(df1diff)
    df1new = pd.concat(df1_parts)


    # --- Handle additional non-hydrogen elements for df2 ---
    df2_parts = [df2common]
    if diff2:
        df2diff = df2NotH.loc[diff2].copy()
        # The starting index for unique elements in df2 should follow the last index used by df1's unique atoms
        # This ensures sequential indexing across both reindexed DFs if they are combined or compared later
        idxMaxFromDf1 = max(df1new.index) if not df1new.empty else -1
        new_diff2_indices = range(idxMaxFromDf1 + 1, idxMaxFromDf1 + 1 + len(diff2))
        df2diff.index = new_diff2_indices
        df2_parts.append(df2diff)
    df2new = pd.concat(df2_parts)

#    print("df2new")
#    print(df2new)
#    print()

    # --- Merge hydrogen atoms back into the reindexed DataFrames ---
    if not df1H.empty:
        # Before concatenating, ensure df1H's indices don't overlap with df1new.
        # This can happen if original df1 had hydrogens with low indices that match reindexed non-H atoms.
        idxMaxCurrentDf1New = max(df1new.index) if not df1new.empty else -1
        # Create new sequential indices for df1H
        new_df1H_indices = range(idxMaxCurrentDf1New + 1, idxMaxCurrentDf1New + 1 + len(df1H))
        df1H_reindexed = df1H.copy()
        df1H_reindexed.index = new_df1H_indices
        df1new = pd.concat([df1new, df1H_reindexed])
    
    if not df2H.empty:
        # Similarly, ensure df2H's indices don't overlap with df2new.
        idxMaxCurrentDf2New = max(df2new.index) if not df2new.empty else -1
        # Create new sequential indices for df2H
        new_df2H_indices = range(idxMaxCurrentDf2New + 1, idxMaxCurrentDf2New + 1 + len(df2H))
        df2H_reindexed = df2H.copy()
        df2H_reindexed.index = new_df2H_indices
        df2new = pd.concat([df2new, df2H_reindexed])

#    print("df2new")
#    print(df2new)
#    print()
    
    # Sort by index to ensure consistent order
    df1new = df1new.sort_index().reset_index(drop=True)
    df2new = df2new.sort_index().reset_index(drop=True)

#    print("df2new")
#    print(df2new)
#    print()

    # Reset ATOM_ID to match new row order (1-based)
    df1new['ATOM_ID'] = range(1, len(df1new) + 1)
    df2new['ATOM_ID'] = range(1, len(df2new) + 1)

    return df1new, df2new

def pad_to_match(df1: pd.DataFrame, df2: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Ensure both DataFrames have the same number of rows by padding the shorter one with empty rows.

    Args:
        df1 (pd.DataFrame): First DataFrame (e.g., reference).
        df2 (pd.DataFrame): Second DataFrame (e.g., referee).

    Returns:
        tuple[pd.DataFrame, pd.DataFrame]: A tuple of the two DataFrames, padded to have the same length.
    """
    if len(df1) != len(df2):
        dfShorter = df2 if len(df1) > len(df2) else df1
        noRowsMissing = abs(len(df1) - len(df2))
        idxPad = range(dfShorter.index[-1] + 1, dfShorter.index[-1] + 1 + noRowsMissing)
        dfPad = pd.DataFrame(index=idxPad, columns=dfShorter.columns)
        if len(df1) > len(df2):
            df2pda = pd.concat([dfShorter, dfPad])
            return df1, df2pda
        else:
            df1pad = pd.concat([dfShorter, dfPad])
            return df1pad, df2
    return df1, df2

#############################################################################################################

def reindex(reference: str, referee: str, outDir: str, visualise=False, logger=None) -> Optional[str]:
    """
    Loads reference molecule (index maintained) and referee molecule (reindexed to match reference).
    Reindexes referee to match reference using RDKit's MCS, logging specific failures and returning None if unsuccessful.

    Args:
        reference (str): Absolute path to a reference molecule structure file.
        referee (str): Absolute path to a referee molecule structure file.
        outDir (str): Directory to save reindexed files.
        reindex_output_name: Identifier for output file naming.
        visualise (bool): Whether to visualize (not implemented).
        print (bool): Whether to print (not implemented).
        logger: Logger instance for error reporting.

    Returns:
        Optional[str]: Path to reindexed PDB file, or None if reindexing fails.

    Raises:
        Exceptions are caught internally and logged; None is returned on failure.
    """
    os.makedirs(outDir, exist_ok=True)

    # Load structure files into dataframe, molecule, name
    try:
        (df1, mol1, name1) = load_molecule(reference)
    except Exception as e:
        if logger:
            logger.error(f"Failed to load reference molecule {reference}: {str(e)}")
        return None

    try:
        (df2, mol2, name2) = load_molecule(referee)
    except Exception as e:
        if logger:
            logger.error(f"Failed to load referee molecule {referee}: {str(e)}")
        return None

    # Get MCS, and indices of atoms matching it and unique to each molecule
    try:
        mcs_mol = get_maximum_common_substructure(mol1, mol2)
    except Exception as e:  # Note: Original code raises SubstructureNotFound, assuming it's an Exception subclass
        if logger:
            logger.error(f"Could not find common substructure between {name1} and {name2}: {str(e)}")
        return None

    try:
        match1, match2 = mol1.GetSubstructMatch(mcs_mol), mol2.GetSubstructMatch(mcs_mol)
        diff1, diff2 = get_atom_difference(mol1, match1), get_atom_difference(mol2, match2)
    except Exception as e:
        if logger:
            logger.error(f"Failed to compute substructure matches or differences for {name2}: {str(e)}")
        return None
    
#    print("match1")
#    print(match1)
#    print()

#    print("match2")
#    print(match2)
#    print()

#    print("diff1")
#    print(diff1)
#    print()

#    print("diff2")
#    print(diff2)
#    print()

    # Reindex the referee dataframe to match atom indices of reference
    try:
        df1new, df2new = get_reindexed_dataframes(df1, df2, match1, match2, diff1, diff2)
    except Exception as e:
        if logger:
            logger.error(f"Failed to reindex dataframes for {name2}: {str(e)}")
        #    print(e)
        return None

    # Save to xyz and pdb files
    referee_reindexed = os.path.join(outDir, f"{name2}_reidx")
    try:
        #df2xyz(df2new, f"{referee_reindexed}.xyz")
        pdbUtils.df2pdb(df2new, f"{referee_reindexed}.pdb")
    except Exception as e:
        if logger:
            logger.error(f"Failed to save reindexed files for {name2}: {str(e)}")
        #    print(e)
        return None

    if visualise:
        save_image_difference(mol1, match1, mol2, match2, f"{outDir}/img_preindexed")

        # Reload reindexed files to visualise the results of reindexing and optimise geometry (only mol needed)
        (_, mol1reidx, _),  (_, mol2reidx, _) = load_molecule(f"{os.path.join(outDir, name1+'_reidx.xyz')}"), load_molecule(f"{os.path.join(outDir, name2+'_reidx.xyz')}")

        # Optimise geometry
        isConv1 = RDforceFields.UFFOptimizeMolecule(mol1reidx)
        if isConv1 != 0:
            raise StructureNotOptimised("Reindexed reference could not have been optimised.")
        RDchem.MolToXYZFile(mol1reidx, f"{os.path.join(outDir, name1+'_reidx_opt.xyz')}")
        isConv2 = RDforceFields.UFFOptimizeMolecule(mol2reidx)
        if isConv2 != 0:
            raise StructureNotOptimised("Reindexed referee could not have been optimised.")
        RDchem.MolToXYZFile(mol2reidx, f"{os.path.join(outDir, name2+'_reidx_opt.xyz')}")
        
        # Get reindexed MCS, and indices of atoms matching it and unique to each molecule
        mcs_molreidx = get_maximum_common_substructure(mol1reidx, mol2reidx)
        match1reidx, match2reidx = mol1reidx.GetSubstructMatch(mcs_molreidx), mol2reidx.GetSubstructMatch(mcs_molreidx)
        # Save an image visualising reindexed atom indices and which atoms are not common to the structures
        save_image_difference(mol1reidx, match1reidx, mol2reidx, match2reidx, f"{outDir}/img_reindexed")
        
        # Provide hint on charge and multiplicty with RDkit
        logger.debug(f"Reference molecule: {name1}")
        logger.debug(f"    Number of unpaired electrons: {Descriptors.NumRadicalElectrons(mol1reidx)}")
        logger.debug(f"    Formal chagre: {RDchem.rdmolops.GetFormalCharge(mol1reidx)}")
        logger.debug(f"Referee molecule: {name2}")
        logger.debug(f"    Number of unpaired electrons: {Descriptors.NumRadicalElectrons(mol2reidx)}")
        logger.debug(f"    Formal chagre: {RDchem.rdmolops.GetFormalCharge(mol2reidx)}")
        
    return f"{referee_reindexed}.pdb"

# Entry point for command-line execution
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a reindexed structure file of the referee to match atom indices and labels to reference input")
    parser.add_argument("--reference", type=str, help="Path to structure file of a molecule whose atom indices are matched")
    parser.add_argument("--referee", type=str, help="Path to structure file of a molecule to be reordered")
    parser.add_argument("--outDir", type=str, help="Path to directory to which output will be saved")
    args = parser.parse_args()
    reference = args.reference
    referee = args.referee
    outDir = args.outDir

    try:
        reindex(reference, referee, outDir, None)
    except (ValueError, TypeError, RuntimeError) as e:
        print(f"Error: {e}")