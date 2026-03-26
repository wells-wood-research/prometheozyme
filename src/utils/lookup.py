import string

class IDLookup:
    def __init__(self, alphabet = "".join(c for c in string.ascii_letters + string.digits if c != "x")):
        # Default: base62
        self.alphabet = alphabet or (string.ascii_lowercase + string.ascii_uppercase + string.digits)
        self.base = len(self.alphabet)

        self.forward = {}   # long_id -> short_id
        self.reverse = {}   # short_id -> long_id
        self.counter = 0

    def _encode(self, num: int) -> str:
        """Convert integer to base-N string."""
        if num == 0:
            return self.alphabet[0]

        chars = []
        while num > 0:
            num, rem = divmod(num, self.base)
            chars.append(self.alphabet[rem])
        return "".join(reversed(chars))

    def get(self, long_id: str) -> str:
        """Get or assign short ID."""
        if long_id in self.forward:
            return self.forward[long_id]

        short_id = self._encode(self.counter)
        self.counter += 1

        self.forward[long_id] = short_id
        self.reverse[short_id] = long_id

        return short_id
    
    def get_site_token(self, site):
        mol_token = self.get(site.molecule_id)
        atom_token = self.get(site.atom_id)
        return f"{mol_token}{atom_token}"

    def resolve(self, short_id: str) -> str:
        """Reverse lookup."""
        return self.reverse[short_id]

    def __contains__(self, long_id: str):
        return long_id in self.forward
    
def build_flavour_map(row, idx_to_flavour):
    mapping = {}
    for i, site in enumerate(row):
        if site is None:
            continue
        flavour_id = idx_to_flavour[i]
        mapping[flavour_id] = (site.molecule_id, site.atom_id)
    return mapping