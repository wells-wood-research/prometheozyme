import os
import uuid
import copy
import pdbUtils
import sys

class Ingredient:
    def __init__(self, path, charge, multiplicity, roles=None, constraints=[], role_title=None, name=None, conformations=None):
        self.path = path
        self.charge = charge
        self.multiplicity = multiplicity
        pdb = pdbUtils.pdb2df(path)
        pdb["ROLE"] = [[] for _ in range(len(pdb))]
        for role_name, atom_names in roles.items():
            mask = pdb["ATOM_NAME"].isin(atom_names)
            pdb.loc[mask, "ROLE"] = pdb.loc[mask, "ROLE"].apply(lambda lst: lst + [role_name])
        self.indices = [(i, row["ATOM_NAME"], row["ROLE"]) for i, row in pdb.iterrows()]
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

def get_constraints_idx(guest, host, constraints):
    for constraint in constraints:
        if constraint.guestIdx in guest.roles.keys():
            print(f"Yes for guest:")
            print(f"{constraint.guestIdx}, {guest.roles[constraint.guestIdx]}")
        if constraint.hostIdx in host.roles.keys():
            print(f"Yes for host:")
            print(f"{constraint.hostIdx}, {host.roles[constraint.hostIdx]}")

    sys.exit(1)
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

def update_constraints(guest, host, constraint, role_title):
    guest_copy = copy.deepcopy(guest)
    guest_copy.id = str(uuid.uuid4())  # Assign a new unique ID
    return guest_copy.update_constraints(host, constraint, role_title)

def update_guest_constraints(roles):
    updated_guests = []
    for role in roles:
        for guest in role.guests:
            role_title = role.title
            constraints = role.constraints
            if constraints:
                for constraint in constraints:
                    updated_guests.append(update_constraints(guest, role.host, constraint, role_title))
            else:
                guest_copy = copy.deepcopy(guest)
                guest_copy.id = str(uuid.uuid4()) 
                updated_guests.append(guest_copy.update_role_title(role_title))
    return updated_guests

def reduce_guests(guests):
    unique_guests = {}
    for guest in guests:
        # Convert constraints to a hashable form by turning lists into tuples
        hashable_constraints = []
        for constraint in guest.constraints:
            guestIdx, guestType, hostIdx, hostType, val = constraint
            # Convert lists to tuples, leave other types unchanged
            hashable_guestIdx = tuple(guestIdx) if isinstance(guestIdx, list) else guestIdx
            hashable_hostIdx = tuple(hostIdx) if isinstance(hostIdx, list) else hostIdx
            hashable_constraints.append((hashable_guestIdx, guestType, hashable_hostIdx, hostType, val))
        key = (guest.path, tuple(hashable_constraints))
        if key not in unique_guests:
            unique_guests[key] = guest
    return list(unique_guests.values())

def print_reduced(updated_guests, unique_guests, logger=None):
    """Print the differences between updated guests and unique guests."""
    missing = set(updated_guests) - set(unique_guests)
    
    if not missing:
        logger.debug("After reducing, no guests were dropped.")
        return
    
    logger.debug("After reducing, the following guests were dropped:")
    for guest in missing:
        logger.debug("Path: %s", guest.path)
        logger.debug("Constraints: %s", guest.constraints)
    logger.debug(f"Total dropped duplicates: {len(missing)}\n")