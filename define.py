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
    def __init__(self, path, eopt, einter, charge, multiplicity, roles=None, constraints=[], role_title=None, name=None, conformations=None, df=None):
        self.path = path
        self.eopt = eopt
        self.einter = einter
        self.charge = charge
        self.multiplicity = multiplicity
        if path.endswith(".pdb") and df is None:
            df = pdbUtils.pdb2df(path)
            df["ROLE"] = [[] for _ in range(len(df))] # TODO rename to flavour
            for role_name, atom_names in roles.items():
                mask = df["ATOM_NAME"].isin(atom_names)
                df.loc[mask, "ROLE"] = df.loc[mask, "ROLE"].apply(lambda lst: lst + [role_name])
            df["ING"] = name
            df["DISH"] = "init"
            df = df[col_order].astype(col_types)
        self.df = df
        self.n_atoms = len(self.df)
        self.constraints = constraints
        self.name = name or os.path.splitext(os.path.basename(path))[0]
        self.id = str(uuid.uuid4())

    def rewrite_xyz(self):
        """Rewrite the XYZ file with updated charge and multiplicity."""
        with open(self.path, 'r') as file:
            lines = file.readlines()
        if len(lines) < 2:
            raise ValueError(f"XYZ file {self.path} does not contain at least two lines.")
        lines[1] = f"{self.charge} {self.multiplicity}\n"
        with open(self.path, 'w') as file:
            file.writelines(lines)

    def get_format(self):
        """Return file format based on extension."""
        return 'pdbfile' if os.path.splitext(self.path)[1] == '.pdb' else 'xyzfile'
    
    def update_constraints(self, host, constraint, role_title):
        self.role_title = role_title
        val = constraint.val
        guestIdx = constraint.guestIdx
        guestType = constraint.guestType
        hostIdx = constraint.hostIdx
        hostType = constraint.hostType

        # Helper function to get all atom indices from an XYZ file
        # Edit: cannot keep all atoms, as non-polar hydrogens are removed by gnina
        def get_all_atom_indices(ingredient):
            try:
                with open(ingredient.path, 'r') as f:
                    atom_count = int(f.readline().strip())
                    f.readline()  # Skip the comment line
                    non_hydrogen_indices = []
                    for i in range(atom_count):
                        line = f.readline().strip()
                        if line:
                            atom_type = line.split()[0]
                            if atom_type.upper() != 'H':
                                non_hydrogen_indices.append(i)
                    return non_hydrogen_indices
            except FileNotFoundError:
                print(f"XYZ file {ingredient.path} not found")
                return []
            except ValueError:
                print(f"Invalid atom count in {ingredient.path}")
                return []

        # Resolve guestIdx
        if guestIdx == "all":
            guestIdx = get_all_atom_indices(self)
        elif isinstance(guestIdx, str):
            guestIdx = getattr(self.indices, guestIdx, guestIdx)
        # Resolve hostIdx
        if hostIdx is None:
            hostIdx = get_all_atom_indices(host)
        elif isinstance(hostIdx, str):
            hostIdx = getattr(host.indices, hostIdx, hostIdx)

        # Check for unresolved symbolic indices
        if isinstance(guestIdx, str) or isinstance(hostIdx, str):
            print(f"Unresolved symbolic constraint indices: guestIdx={guestIdx}, hostIdx={hostIdx}. Check spelling and presence in indices.")
            return self

        # Ensure guestIdx and hostIdx are lists for consistent processing
        guest_indices = guestIdx if isinstance(guestIdx, list) else [guestIdx]
        host_indices = hostIdx if isinstance(hostIdx, list) else [hostIdx]

        # Create constraints for all combinations of guest and host indices
        keep = (guest_indices, guestType, host_indices, hostType, val)
        if keep not in self.constraints:
            self.constraints.append(keep)
        return self
    
    def update_role_title(self, role_title):
        """Update the role title of the ingredient."""
        self.role_title = role_title
        return self
    
    def update_conformations(self, conformations):
        """Update the role title of the ingredient."""
        self.conformations = conformations
        return self

class Role:
    def __init__(self, title, priority, guests, host, constraints=None):
        self.title = title
        self.priority = priority
        self.guests = guests
        self.host = host
        self.constraints = [] if constraints is None else (constraints if isinstance(constraints, list) else [constraints])

class Constraint:
    def __init__(self, guestIdx, hostIdx, val, guestType="any", hostType="any", force=100):
        self.guestIdx = guestIdx
        self.guestType = guestType # TODO remove? Allow com or iter?
        self.hostIdx = hostIdx
        self.hostType = hostType
        self.val = val
        self.force = force