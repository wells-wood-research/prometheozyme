import argparse

def read_xyz(file_path):
    """Read an XYZ file and return a list of (atom_count, comment, coordinates)."""
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    structures = []
    i = 0
    while i < len(lines):
        try:
            atom_count = int(lines[i].strip())
            comment = lines[i + 1].strip()
            parts = [line.strip().split() for line in lines[i + 2:i + 2 + atom_count]]
            coordinates = [f"{part[0]}{float(part[1]):27.17f}{float(part[2]):27.17f}{float(part[3]):27.17f}\n" for part in parts]
            structures.append((atom_count, comment, coordinates))
            i += atom_count + 2
        except (ValueError, IndexError):
            break  # Stop if no more valid structures
    return structures

def merge_xyz(receptor_path, ligand_path, output_path):
    """Merge receptor XYZ with each structure in a multi-XYZ ligand file."""
    # Read receptor (single structure)
    receptor_structures = read_xyz(receptor_path)
    if not receptor_structures:
        raise ValueError("Receptor file is empty or invalid")
    receptor_count, receptor_comment, receptor_coords = receptor_structures[0]
    
    # Read ligand (potentially multi-XYZ)
    ligand_structures = read_xyz(ligand_path)
    if not ligand_structures:
        raise ValueError("Ligand file is empty or invalid")
    
    # Write merged multi-XYZ file
    with open(output_path, 'w') as f:
        for ligand_count, ligand_comment, ligand_coords in ligand_structures:
            # Calculate total atoms for this structure
            total_atoms = receptor_count + ligand_count
            # Create comment for merged structure
            merged_comment = f"Merged from {receptor_path} and {ligand_path} (structure)"
            # Write structure block
            f.write(f"{total_atoms}\n")
            f.write(f"{merged_comment}\n")
            f.writelines(receptor_coords)
            f.writelines(ligand_coords)
            f.write("\n")  # Separator between structures

def main(receptor, ligand, output):
    merge_xyz(receptor, ligand, output)

    print(f"Merged XYZ file written to {output}")

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Merge two XYZ files into one.")
    parser.add_argument("--receptor", required=False, help="Path to receptor XYZ file")
    parser.add_argument("--ligand", required=False, help="Path to ligand XYZ file")
    parser.add_argument("--output", required=False, help="Path to output merged XYZ file")
    
    # Parse arguments
    args = parser.parse_args()
    receptor, ligand, output = args.receptor, args.ligand, args.output
    if not all([receptor, ligand, output]):
        receptor = "/home/mchrnwsk/theozymes/ingredients/sub+H2./0_SUB_opt_unrest/SUB_opt_unrest.xyz"
        ligand = "/home/mchrnwsk/theozymes/docking/covalent/out.xyz"
        output = "/home/mchrnwsk/theozymes/docking/covalent/out_merged.xyz"        
    
    # Merge the files
    main(receptor, ligand, output)