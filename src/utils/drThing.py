import os
import uuid
import logging
import pandas as pd
from typing import List, Dict, Optional

from utils.drStructure import pdb2df

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
## CUSTOM CLASSES
########################

col_order = [
    "ATOM", "ATOM_ID", "ATOM_NAME", "RES_NAME", "CHAIN_ID", "RES_ID",
    "X", "Y", "Z", "OCCUPANCY", "BETAFACTOR", "ELEMENT", "FLAVOUR", "ING", "DISH"
]

col_types = {
    "ATOM": str,
    "ATOM_ID": int,
    "ATOM_NAME": str,
    "RES_NAME": str,
    "CHAIN_ID": str,
    "RES_ID": int,
    "X": float,
    "Y": float,
    "Z": float,
    "OCCUPANCY": float,
    "BETAFACTOR": float,
    "ELEMENT": str,
    "FLAVOUR": str,
    "ING": str,
    "DISH": str
}

class Ingredient:
    def __init__(
        self,
        pathPDB: str,
        pathXYZ: str,
        name: Optional[str],
        eopt: float = 0.0,
        einter: float = 0.0,
        charge: int = 0,
        multiplicity: int = 1,
        flavours: Optional[Dict[str, List[str]]] = None,
        df: Optional[pd.DataFrame] = None,
    ) -> None:
        self.pathPDB: str = pathPDB
        self.pathXYZ: str = pathXYZ
        self.name: str = name or os.path.splitext(os.path.basename(pathPDB))[0]
        self.eopt: float = eopt
        self.einter: float = einter
        self.charge: int = int(charge)
        self.multiplicity: int = int(multiplicity)

        # DataFrame setup
        if df is None:
            df = pdb2df(self.pathPDB)
        if df is None:
            raise ValueError(f"Failed to load PDB data from {self.pathPDB}")

        # Fill in missing columns
        if "FLAVOUR" not in df.columns:
            df["FLAVOUR"] = [[] for _ in range(len(df))]
            if flavours:
                for role_name, atom_names in flavours.items():
                    mask = df["ATOM_NAME"].isin(atom_names)
                    df.loc[mask, "FLAVOUR"] = df.loc[mask, "FLAVOUR"].apply(
                        lambda lst: lst + [role_name]
                    )

        if "ING" not in df.columns:
            df["ING"] = self.name
        if "DISH" not in df.columns:
            df["DISH"] = "init"

        # col_order and col_types are assumed global
        df = df[col_order].astype(col_types)

        self.df: pd.DataFrame = df
        self.n_atoms: int = len(self.df)
        self.id: str = str(uuid.uuid4())


class Course:
    def __init__(
        self,
        name: str,
        host: Ingredient,
        guests: List[Ingredient],
        restraints: Optional[List["Restraint"]] = None,
    ) -> None:
        self.name: str = name
        self.host: Ingredient = host
        self.guests: List[Ingredient] = guests
        self.restraints: List[Restraint] = restraints


class Restraint:
    def __init__(
        self,
        guestIdx: int,
        hostIdx: int,
        val: float,
        tol: float,
        force: float = 100.0,
    ) -> None:
        self.guestIdx: int = guestIdx
        self.hostIdx: int = hostIdx
        self.val: float = val
        self.tol: float = tol
        self.force: float = force
