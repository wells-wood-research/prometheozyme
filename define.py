import os
import uuid
import copy
import pdbUtils
import sys

col_order = [
    "ATOM", "ATOM_ID", "ATOM_NAME", "RES_NAME", "CHAIN_ID", "RES_ID",
    "X", "Y", "Z", "OCCUPANCY", "BETAFACTOR", "ELEMENT", "ROLE", "ING", "DISH"
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
    "ROLE": str,
    "ING": str,
    "DISH": str
}

class Ingredient:
    def __init__(self, name, path, charge=0, multiplicity=1, eopt=0, einter=0,  flavours=None, restraints=[], df=None):
        self.path = path
        self.eopt = eopt
        self.einter = einter
        self.charge = int(charge)
        self.multiplicity = int(multiplicity)
        # TODO consistency with PDB paths
        if path.endswith(".pdb") and df is None:
            df = pdbUtils.pdb2df(path)
            df["ROLE"] = [[] for _ in range(len(df))] # TODO rename to flavour
            for role_name, atom_names in flavours.items():
                mask = df["ATOM_NAME"].isin(atom_names)
                df.loc[mask, "ROLE"] = df.loc[mask, "ROLE"].apply(lambda lst: lst + [role_name])
            df["ING"] = name
            df["DISH"] = "init"
            df = df[col_order].astype(col_types)
        self.df = df
        self.n_atoms = len(self.df)
        self.restraints = restraints
        self.name = name or os.path.splitext(os.path.basename(path))[0]
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