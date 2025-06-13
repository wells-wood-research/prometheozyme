from pdbUtils import pdbUtils
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

from utils.errors import XYZFileFormatError, SubstructureNotFound, StructureNotOptimised

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
        RDchem.rdDetermineBonds.DetermineBonds(mol,charge=0) # TODO Must allow user to input their charge
    mol.SetProp("name", name)
    return df, mol, name