from pymol import cmd
import os

def protonate(path=None, unprotonatedIdx=None):
    # Load your molecule (replace 'molecule.pdb' with your file)
    cmd.load(path, 'mol')

    # Add hydrogens to the molecule
    cmd.h_add('mol')

    # Define the list of atom indices (example: 10, 15, 20) to remove
    if unprotonatedIdx:
        # Convert orca indices to a PyMOL selection string
        pymol_indices = [x+1 for x in unprotonatedIdx]
        index_selection = '+'.join(str(i) for i in pymol_indices)
        # Note: Use 'neighbor' to select hydrogens bonded to the target atoms
        cmd.remove(f'elem H and neighbor (index {index_selection})')

    # Save the modified structure with hydrogens (optional)
    outfile = f"{os.path.splitext(path)[0]}_Hs.pdb"
    cmd.save(outfile, 'mol', state=0, format='pdb', quiet=0)

    return outfile