import numpy as np
import argparse
from typing import Tuple, List
import subprocess
import os

# Van der Waals radii in Angstroms for common elements
VDW_RADII = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80, 'P': 1.80,
    'F': 1.47, 'Cl': 1.75, 'Br': 1.85, 'I': 1.98
}

def read_xyz_file(file_path: str) -> Tuple[List[str], np.ndarray]:
    """
    Read an XYZ file and return atom types and their coordinates.
    
    Args:
        file_path: Path to the XYZ file
        
    Returns:
        Tuple containing list of atom types and numpy array of coordinates
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    atom_types = []
    coords = []
    for line in lines[2:]:  # Skip atom count and comment lines
        parts = line.strip().split()
        if len(parts) >= 4:
            atom_types.append(parts[0])
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    
    return atom_types, np.array(coords)

def calculate_center_and_extent(coords: np.ndarray, atom_types: List[str]) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the center and extent of a molecule using its spatial extent with VdW radii.
    
    Args:
        coords: Array of atom coordinates (n x 3)
        atom_types: List of atom types
        
    Returns:
        Tuple of (center_x, center_y, center_z, size_x, size_y, size_z)
    """
    vdw_radii = [VDW_RADII.get(t, 1.5) for t in atom_types]  # Default 1.5Å
    max_radii = max(vdw_radii)
    min_coords = np.min(coords, axis=0) - max_radii
    max_coords = np.max(coords, axis=0) + max_radii
    
    cx = (max_coords[0] + min_coords[0]) / 2
    cy = (max_coords[1] + min_coords[1]) / 2
    cz = (max_coords[2] + min_coords[2]) / 2
    sx = max_coords[0] - min_coords[0] / 2
    sy = max_coords[1] - min_coords[1] / 2
    sz = max_coords[2] - min_coords[2] /2 
    
    return cx, cy, cz, sx, sy, sz

def calculate_docking_box(receptor_path: str, ligand_path: str, padding: float = 10.0) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate docking box parameters based on receptor and ligand XYZ files.
    
    Args:
        receptor_path: Path to receptor XYZ file
        ligand_path: Path to ligand XYZ file
        padding: Extra space around receptor in Angstroms (default: 10.0)
        
    Returns:
        Tuple of (center_x, center_y, center_z, size_x, size_y, size_z)
    """
    # Read receptor and ligand files
    receptor_types, receptor_coords = read_xyz_file(receptor_path)
    ligand_types, ligand_coords = read_xyz_file(ligand_path)
    
    # Calculate receptor center and extent
    center_x, center_y, center_z, receptor_size_x, receptor_size_y, receptor_size_z = calculate_center_and_extent(receptor_coords, receptor_types)
    
    # Calculate ligand extent
    _, _, _, ligand_size_x, ligand_size_y, ligand_size_z = calculate_center_and_extent(ligand_coords, ligand_types)
    
    # Ligand's maximum dimension
    ligand_max_dim = max(ligand_size_x, ligand_size_y, ligand_size_z)
    
    # Box size: receptor extent + ligand max dimension + padding
    size_x = receptor_size_x + ligand_max_dim + 2 * padding
    size_y = receptor_size_y + ligand_max_dim + 2 * padding
    size_z = receptor_size_z + ligand_max_dim + 2 * padding
    
    return center_x, center_y, center_z, size_x, size_y, size_z

def main(receptor, ligand, padding):
    # Calculate docking box
    center_x, center_y, center_z, size_x, size_y, size_z = calculate_docking_box(
        receptor, ligand, padding
    )
    
    # Print results
    print(f"Docking Box Parameters:")
    print(f"  center_x: {center_x:.3f} Å")
    print(f"  center_y: {center_y:.3f} Å")
    print(f"  center_z: {center_z:.3f} Å")
    print(f"  size_x: {size_x:.3f} Å")
    print(f"  size_y: {size_y:.3f} Å")
    print(f"  size_z: {size_z:.3f} Å")

def dock(host, ingredient, outdir, dock_params, redocking=False, logger=None):
    receptor = host.path
    ligand = ingredient.path
    outdir = os.path.join(outdir, ingredient.name)
    os.makedirs(outdir, exist_ok=True)  # Ensure output directory exists
    
    # Calculate docking box
    center_x, center_y, center_z, size_x, size_y, size_z = calculate_docking_box(receptor, ligand, padding=5.0)
    
    # Update dock_params with docking box coordinates
    dock_params = dock_params.copy()  # Avoid modifying the original
    dock_params.update({
        'center_x': center_x,
        'center_y': center_y,
        'center_z': center_z,
        'size_x': size_x,
        'size_y': size_y,
        'size_z': size_z
    })

    # Construct the command as a list
    cmd = [
        "gnina",
        "-r", receptor,
        "-l", ligand,
        "--center_x", str(dock_params['center_x']),
        "--center_y", str(dock_params['center_y']),
        "--center_z", str(dock_params['center_z']),
        "--size_x", str(dock_params['size_x']),
        "--size_y", str(dock_params['size_y']),
        "--size_z", str(dock_params['size_z']),
        "--scoring", dock_params['scoring'],
        "--cnn_scoring", dock_params['cnn_scoring'],
        "--pose_sort_order", dock_params['pose_sort_order'],
        "-o", os.path.join(outdir, "out.xyz"),
        "--atom_terms", os.path.join(outdir, "atom_terms"),
        "--exhaustiveness", str(dock_params['exhaustiveness']),
        "--num_modes", str(dock_params['num_modes'])
    ]

    # Include optional parameters if not None or False
    if dock_params['addH']:
        cmd.append("--addH")
    if dock_params['stripH']:
        cmd.append("--stripH")
    if dock_params.get("no_gpu", False):
       cmd.append("--no_gpu")
    if dock_params.get("min_rmsd_filter", None):
        cmd.extend(["--min_rmsd_filter", str(dock_params['min_rmsd_filter'])])
    if dock_params.get("temperature", None):
        cmd.extend(["--temperature", str(dock_params['temperature'])])
    if dock_params.get("num_mc_steps", None):
        cmd.extend(["--num_mc_steps", str(dock_params['num_mc_steps'])])

    with open(os.path.join(outdir, 'scores.txt'), 'w') as outfile:
        subprocess.run(cmd, check=True, stdout=outfile, stderr=subprocess.STDOUT)
    logger.info(f"Docking for {ingredient.name} {'(redocking)' if redocking else ''} completed. Results saved in {outdir}")
    return os.path.join(outdir, "out.xyz")  # Return the path to the docked output file


if __name__ == '__main__':
    """Calculate docking box parameters from receptor and ligand XYZ files."""
    parser = argparse.ArgumentParser(description='Calculate docking box parameters for ligand docking.')
    parser.add_argument('--receptor', type=str, required=False, help='Path to receptor XYZ file')
    parser.add_argument('--ligand', type=str, required=False, help='Path to ligand XYZ file')
    parser.add_argument('--padding', type=float, default=5.0, help='Padding around receptor in Angstroms (default: 5.0)')
    args = parser.parse_args()

    receptor, ligand, padding = args.receptor, args.ligand, args.padding
    if not all([receptor, ligand]):
        receptor = "/home/mchrnwsk/theozymes/ingredients/sub+H2./0_SUB_opt_unrest/SUB_opt_unrest.xyz"
        ligand = "/home/mchrnwsk/theozymes/docking/covalent/out.xyz"
    
    main(receptor, ligand, padding)