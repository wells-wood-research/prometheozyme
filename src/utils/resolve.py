from utils import structure, types

class StructureFileResolver:
    def __init__(self, encoder, config, cwd, outdir, rows):
        """
        Parameters
        ----------
        encoder : SignatureEncoder
        config : parsed config (contains ingredients)
        cwd : Path
        outdir : Path
        structure : your structure module (for path helpers)
        """
        self.encoder = encoder
        self.config = config
        self.cwd = cwd
        self.outdir = outdir
        self.structure = structure

        # Precompute library signatures
        self.library = self._build_library_index(rows)

    # -------------------------
    # Public API
    # -------------------------

    def resolve(self, host):
        """
        Resolve host motif → structure file(s)
        """

        sites = [site for _, site in host]
        mol_ids = {s.molecule_id for s in sites}

        # Case 1: single molecule → use library
        if len(mol_ids) == 1:
            mol_id = next(iter(mol_ids))

            return [
                structure.get_molec_xyz_path(
                    self.cwd,
                    self.config.ingredients[mol_id].filepath
                )
            ]

        # Case 2: multi-molecule → assembly
        host_signature = self.encoder.motif_to_signature(host)

        host_dir = structure.prep_assembly_dir(self.outdir, host_signature)

        return [host_dir / "dock.docker.xyz"]

    # -------------------------
    # Internal
    # -------------------------

    def _build_library_index(self, rows):
        """
        Build: signature → molecule_id
        using actual recipe rows
        """
        index = {}

        for row in rows:
            # Convert row (tuple of Sites) → motif
            motif = frozenset((i, site) for i, site in enumerate(row))

            signature = self.encoder.motif_to_signature(motif)

            # Extract molecule identity
            mol_ids = {site.molecule_id for site in row}

            if len(mol_ids) == 1:
                mol_id = next(iter(mol_ids))
                if len(mol_ids) == 1:
                    # ensure it's a FULL molecule (not partial motif)
                    expected_cols = {
                        i for i, site in enumerate(row)
                        if site.molecule_id == mol_id
                    }

                    motif_cols = {i for i, _ in motif}

                    if motif_cols == expected_cols:
                        index[signature] = mol_id

        return index

    def _make_site(self, mol_id, atom_id):
        """
        Factory for Site objects (adapt to your types module)
        """
        return types.Site(molecule_id=mol_id, atom_id=atom_id)