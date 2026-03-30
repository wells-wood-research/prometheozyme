class SignatureEncoder:
    def __init__(self, lookup, n_cols, empty_token="x", sep="-"):
        """
        Parameters
        ----------
        lookup : object
            Must provide: get_site_token(site)
        n_cols : int
            Total number of columns in the motif/row
        empty_token : str
            Token used for empty positions
        sep : str
            Separator for signature string
        """
        self.lookup = lookup
        self.n_cols = n_cols
        self.empty_token = empty_token
        self.sep = sep

    # -------------------------
    # Core transformations
    # -------------------------

    def motif_to_row(self, motif):
        """
        Convert motif (sparse) → row (dense list)

        motif: iterable of (col, Site)
        """
        row = [None] * self.n_cols
        for c, site in motif:
            row[c] = site
        return row

    def row_to_signature(self, row):
        """
        Convert row → signature string
        """
        tokens = [
            self.empty_token if site is None
            else self.lookup.get_site_token(site)
            for site in row
        ]
        return self.sep.join(tokens)

    def motif_to_signature(self, motif):
        """
        Direct motif → signature
        """
        row = self.motif_to_row(motif)
        return self.row_to_signature(row)