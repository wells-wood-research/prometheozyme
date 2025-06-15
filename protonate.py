from pymol import cmd
import os

def protonate_all(path=None):
    # Clear all objects in PyMOL
    cmd.delete('all')

    cmd.set("retain_order", "on")

    # Load your molecule
    cmd.load(path, 'mol')

    # Add hydrogens to the molecule
    cmd.h_add('mol')

    # Save the fully protonated structure
    outfile = f"{os.path.splitext(path)[0]}_prot.pdb"
   
    cmd.save(outfile, 'mol', state=0, format='pdb', quiet=0)
    atom_count = cmd.count_atoms('mol')

    # Delete the molecule to clean up
    cmd.delete('mol')

    return outfile, atom_count

def deprotonate_selected(path=None, unprotonatedIdx=None):
    # Clear all objects in PyMOL
    cmd.delete('all')

    # Load your molecule
    cmd.load(path, 'mol')

    # Remove hydrogens from specified atoms
    if unprotonatedIdx:
        # Convert 0-based indices to 1-based PyMOL indices
        pymol_indices = [x + 1 for x in unprotonatedIdx]
        index_selection = '+'.join(str(i) for i in pymol_indices)
        cmd.remove(f'elem H and neighbor (index {index_selection})')

    # Save the deprotonated structure
    outfile = f"{os.path.splitext(path)[0]}_deprot.pdb"
    cmd.set("retain_order", "on")
    cmd.save(outfile, 'mol', state=0, format='pdb', quiet=0)
    atom_count = cmd.count_atoms('mol')

    # Delete the molecule to clean up
    cmd.delete('mol')

    return outfile, atom_count