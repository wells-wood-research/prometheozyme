import os
import uuid
from utils import pdb2df

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
    def __init__(self, pathPDB, pathXYZ, name, eopt=0, einter=0, charge=0, multiplicity=1, flavours=None, df=None):
        self.pathPDB = pathPDB
        self.pathXYZ = pathXYZ
        self.name = name or os.path.splitext(os.path.basename(pathPDB))[0]
        self.eopt = eopt
        self.einter = einter
        self.charge = int(charge)
        self.multiplicity = int(multiplicity)
        # pathPDB and df could contain different information so don't overwrite df if it's passed in
        if df is None:
            df = pdb2df(self.pathPDB)
        # Fill in missing df columns
        if "FLAVOUR" not in df.columns:
            df["FLAVOUR"] = [[] for _ in range(len(df))]
            for role_name, atom_names in flavours.items():
                mask = df["ATOM_NAME"].isin(atom_names)
                df.loc[mask, "FLAVOUR"] = df.loc[mask, "FLAVOUR"].apply(lambda lst: lst + [role_name])
        if "ING" not in df.columns:
            df["ING"] = self.name
        if "DISH" not in df.columns:
            df["DISH"] = "init"
        df = df[col_order].astype(col_types)
        self.df = df
        self.n_atoms = len(self.df)
        self.id = str(uuid.uuid4())

class Course:
    def __init__(self, name, host, guests, restraints=None):
        self.name = name
        self.host = host
        self.guests = guests
        self.restraints = [] if restraints is None else (restraints if isinstance(restraints, list) else [restraints])

class Restraint:
    def __init__(self, guestIdx, hostIdx, val, force=100):
        self.guestIdx = guestIdx
        self.hostIdx = hostIdx
        self.val = val
        self.force = force