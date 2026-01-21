import ast
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
    "ING": str,
    "DISH": str
}

class Selection:    
    def __init__(
        self,
        parent: str,
        idx: str
    ) -> None:
        self.parent: str = parent
        self.idx: str = idx
    def __str__(self):
        return f"Selection(parent={self.parent}, idx={self.idx})"

    @classmethod
    def from_dict(cls, d):
        return cls(d["parent"], d["idx"])
    
def _abs_index(selection: Selection, host_n_atoms: int) -> int:
    """Return absolute index for a Selection in merged host+guest system."""
    if selection.parent == "host":
        return int(selection.idx)
    elif selection.parent == "guest":
        return int(selection.idx) + int(host_n_atoms)
    else:
        raise ValueError(f"Unknown selection parent: {selection.parent}")

class Params:
    def __init__(
        self,
        val: float,
        uptol: Optional[float] = None,
        downtol: Optional[float] = None,
        force: float = 100.0,
    ) -> None:
        self.val: float = val
        self.uptol: Optional[float] = uptol
        self.downtol: Optional[float] = downtol
        self.force: float = force
    def __str__(self):
        return f"Params(val={self.val}, uptol={self.uptol}, downtol={self.downtol}, force={self.force})"
    @classmethod
    def from_dict(cls, d):
        return cls(d["val"], d.get("uptol", None), d.get("downtol", None), d.get("force", 100))

class Restraint:    
    def __init__(
        self,
        property: str,
        sele: List[Selection],
        params: Params,
        step: str
    ) -> None:
        self.property: str = property
        self.sele: List[Selection] = [Selection.from_dict(a) for a in sele]
        self.params: Params = Params.from_dict(params)
        self.step: str = step
    def __str__(self):
        sele_str = ", ".join(str(s) for s in self.sele)
        return (
            f"Restraint(property={self.property}, "
            f"atoms=[{sele_str}], "
            f"params={self.params})"
        )

def ing_dish_signature(df):
    pairs = []
    last = None

    for ing, dish in zip(df["ING"], df["DISH"]):
        curr = (ing, dish)
        if curr != last:
            pairs.append(curr)
            last = curr

    return tuple(pairs)

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
        restraints: List[Restraint] = []
    ) -> None:
        self.pathPDB: str = pathPDB
        self.pathXYZ: str = pathXYZ
        self.name: str = name or os.path.splitext(os.path.basename(pathPDB))[0]
        self.eopt: float = eopt
        self.einter: float = einter
        self.charge: int = int(charge)
        self.multiplicity: int = int(multiplicity)
        self.restraints: List[Restraint] = restraints

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
    def __str__(self):
        flavour_summary = {
            k: len(v) for k, v in (
                {role: self.df[self.df["FLAVOUR"].apply(lambda lst: role in lst)]["ATOM_NAME"].tolist() 
                 for role in self.df["FLAVOUR"].explode().dropna().unique()}
            ).items()
        } if "FLAVOUR" in self.df.columns else {}
        return (
            f"Ingredient: {self.name}\n"
            f"PDB_path: {self.pathPDB}\n"
            f"XYZ_path: {self.pathXYZ}\n"
            f"Charge: {self.charge}, Multiplicity: {self.multiplicity}\n"
            f"Eopt: {self.eopt:.4f}, Einter: {self.einter:.4f}\n"
            f"Number of atoms: {self.n_atoms}\n"
            f"Flavours: {flavour_summary}\n"
        )
    
class Course:
    def __init__(
        self,
        name: str,
        host: Ingredient,
        guests: List[Ingredient],
        restraints: Optional[List[Restraint]] = None,
        orcaSettings: Optional[Dict] = None
    ) -> None:
        self.name: str = name
        self.host: Ingredient = host
        self.guests: List[Ingredient] = guests
        self.restraints: List[Restraint] = restraints
        self.orcaSettings: Optional[Dict] = orcaSettings
    def __str__(self):
        guest_names = [g.name for g in self.guests]
        restraint_summaries = "\n    ".join(str(r) for r in self.restraints) if self.restraints else "None"
        return (
            f"Course: {self.name}\n"
            f"Host: {self.host.name}\n"
            f"Guests: {guest_names}\n"
            f"Restraints:\n    {restraint_summaries}\n"
            f"ORCA Settings: {self.orcaSettings or 'Default'}\n"
        )
