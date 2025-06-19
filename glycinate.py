import pandas as pd
from dock import calculate_docking_box
import pdbUtils
import os

def get_next_chain_id(used_chains):
    """Get the next available chain ID."""
    import string
    available = [c for c in string.ascii_uppercase if c not in used_chains]
    return available[0] if available else 'A'

def add_glycine_to_pdb(input_pdb, output_pdb, glycine_pdb="/home/mchrnwsk/theozymes/ingredients/glycine.pdb"):
    """Add a glycine residue to the PDB file."""
    # Read input PDB and glycine PDB as DataFrames
    structure = pdbUtils.pdb2df(input_pdb)
    glycine = pdbUtils.pdb2df(glycine_pdb)

    # Get maximum atom and residue numbers
    max_atom_num = max(structure["ATOM_ID"])
    max_residue_num = max(structure["RES_ID"])

    # Get next chain ID and residue number
    new_chain_id = get_next_chain_id(set(structure["CHAIN_ID"]))
    new_residue_num = max_residue_num + 1
    new_atom_num = max_atom_num + 1

    # Calculate docking box to determine a safe placement
    center_x, center_y, center_z, size_x, size_y, size_z = calculate_docking_box(input_pdb, input_pdb, padding=10.0)

    # Update glycine DataFrame
    glycine = glycine.copy()  # Avoid modifying the original
    glycine["ATOM_ID"] = range(new_atom_num, new_atom_num + len(glycine))
    glycine["RES_ID"] = new_residue_num
    glycine["CHAIN_ID"] = new_chain_id

    # Translate glycine to a position outside the docking box (at least 10 Å away)
    # Place along x-axis at center_x + size_x/2 + 10 Å
    translation_x = center_x + size_x / 2
    glycine["X"] += translation_x
    glycine["Y"] += center_y
    glycine["Z"] += center_z

    # Concatenate DataFrames
    combined_df = pd.concat([structure, glycine], ignore_index=True)

    # Save to output PDB
    pdbUtils.df2pdb(combined_df, output_pdb)

# Example usage
if __name__ == "__main__":
    input_pdb = "/home/mchrnwsk/theozymes/docking/output_2025-06-18_17-53-23_copy/final_theozymes/arrangement_1.pdb"
    output_pdb = f"{os.path.splitext(input_pdb)[0]}_glycine.pdb"
    add_glycine_to_pdb(input_pdb, output_pdb)
    print(f"Glycine added to {output_pdb}")