import os
import uuid
import copy

class Ingredient:
    def __init__(self, path, charge, multiplicity, indices=None, constraints=[], role_title=None, name=None):
        self.path = path
        self.charge = charge
        self.multiplicity = multiplicity
        self.indices = indices
        self.constraints = constraints
        self.role_title = role_title
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
        def get_all_atom_indices(ingredient):
            try:
                with open(ingredient.path, 'r') as f:
                    atom_count = int(f.readline().strip())
                return list(range(atom_count))  # Indices start at 0
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

class Role:
    def __init__(self, title, priority, guests, host, constraints=None):
        self.title = title
        self.priority = priority
        self.guests = guests
        self.host = host
        self.constraints = [] if constraints is None else (constraints if isinstance(constraints, list) else [constraints])

class Constraint:
    def __init__(self, guestIdx, hostIdx, val, guestType="iter", hostType="iter"):
        self.guestIdx = guestIdx
        self.guestType = guestType
        self.hostIdx = hostIdx
        self.hostType = hostType
        self.val = val

class Indices:
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
            
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