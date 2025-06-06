import numpy as np
import argparse
from typing import List, Tuple, Dict

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
    
    # Skip first two lines (atom count and comment)
    atom_types = []
    coords = []
    for line in lines[2:]:
        parts = line.strip().split()
        if len(parts) >= 4:
            atom_types.append(parts[0])
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    
    return atom_types, np.array(coords)

def create_grid(coords: np.ndarray, vdw_radii: List[float], voxel_size: float = 1.0) -> Tuple[np.ndarray, Tuple[int, int, int]]:
    """
    Create a grid that encompasses the molecule with Van der Waals radii.
    
    Args:
        coords: Array of atom coordinates
        vdw_radii: List of Van der Waals radii for each atom
        voxel_size: Size of each voxel in Angstroms
    
    Returns:
        Tuple of grid dimensions and grid size
    """
    # Find min and max coordinates including VdW radii
    max_radii = max(vdw_radii)
    min_coords = np.min(coords, axis=0) - max_radii
    max_coords = np.max(coords, axis=0) + max_radii

    cx = (max_coords[0] + min_coords[0]) / 2
    sx = (max_coords[0] - min_coords[0]) / 2
    cy = (max_coords[1] + min_coords[1]) / 2
    sy = (max_coords[1] - min_coords[1]) / 2
    cz = (max_coords[2] + min_coords[2]) / 2
    sz = (max_coords[2] - min_coords[2]) / 2    
    
    # Calculate grid size
    # grid_size = np.ceil((max_coords - min_coords) / voxel_size).astype(int)
    # return min_coords, grid_size
    print(cx,cy,cz,sx,sy,sz)
    return cx,cy,cz,sx,sy,sz

def map_molecule_to_grid(atom_types: List[str], coords: np.ndarray, min_coords: np.ndarray, 
                        grid_size: Tuple[int, int, int], voxel_size: float = 1.0) -> np.ndarray:
    """
    Map molecule atoms to a voxel grid, calculating occupancy percentage.
    
    Args:
        atom_types: List of atom types
        coords: Array of atom coordinates
        min_coords: Minimum coordinates of the grid
        grid_size: Size of the grid in each dimension
        voxel_size: Size of each voxel in Angstroms
    
    Returns:
        4D numpy array with voxel occupancy and atom index
    """
    print(f"Grid size: {grid_size}")
    
    # Initialize grid
    try:
        grid = np.zeros(grid_size + (2,), dtype=np.float32)  # [x,y,z,(occupancy, atom_index)]
        print(f"Grid shape: {grid.shape}")  # Debug: Print grid shape
    except TypeError as e:
        print(f"Error initializing grid: {e}")
        raise ValueError(f"Invalid grid_size format: {grid_size}. Must be a tuple of 3 integers.")
    
    for i, (atom, coord, radius) in enumerate(zip(atom_types, coords, 
                                                [VDW_RADII.get(t, 1.5) for t in atom_types])):
        # Convert atom coordinates to grid indices
        grid_coord = ((coord - min_coords) / voxel_size).astype(int)
        
        # Calculate the range of voxels potentially affected
        radius_voxels = int(np.ceil(radius / voxel_size))
        x_range = range(max(0, grid_coord[0] - radius_voxels), 
                       min(grid_size[0], grid_coord[0] + radius_voxels + 1))
        y_range = range(max(0, grid_coord[1] - radius_voxels), 
                       min(grid_size[1], grid_coord[1] + radius_voxels + 1))
        z_range = range(max(0, grid_coord[2] - radius_voxels), 
                       min(grid_size[2], grid_coord[2] + radius_voxels + 1))
        
        # Check each voxel in the range
        for x in x_range:
            for y in y_range:
                for z in z_range:
                    # Calculate voxel center
                    voxel_center = min_coords + np.array([x, y, z]) * voxel_size + voxel_size/2
                    distance = np.linalg.norm(voxel_center - coord)
                    
                    # If voxel center is within VdW radius, assign occupancy
                    if distance <= radius:
                        try:
                            grid[x, y, z, 0] = 1.0
                            grid[x, y, z, 1] = i + 1  # Store 1-based index (line number - 2)
                        except IndexError as e:
                            print(f"IndexError at (x={x}, y={y}, z={z}): {e}")
                            print(f"Grid shape: {grid.shape}")
                            raise
    
    return grid

def main(input_file: str):
    """Main function to process XYZ file and create voxel grid."""
    # Read XYZ file
    atom_types, coords = read_xyz_file(input_file)
    
    # Get Van der Waals radii
    vdw_radii = [VDW_RADII.get(atom, 1.5) for atom in atom_types]  # Default radius 1.5Å
    
    # Create grid
    min_coords, grid_size = create_grid(coords, vdw_radii)
    
    # Map molecule to grid
    grid = map_molecule_to_grid(atom_types, coords, min_coords, grid_size)
    
    # Print results
    print(f"Grid dimensions: {grid_size}")
    print(f"Grid origin: {min_coords}")
    print(f"Number of occupied voxels: {np.sum(grid[:,:,:,0] > 0)}")
    
    # Optionally save the grid (uncomment to use)
    np.save('molecule_grid.npy', grid)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Map molecule from XYZ file to voxel grid.')
    parser.add_argument('input_file', required=False, help='Path to input XYZ file')
    args = parser.parse_args()
    input_file = args.input_file
    if input_file is None:
        input_file = "/home/mchrnwsk/theozymes/ingredients/lmf/0_LMF_opt_unrest/LMF_opt_unrest.xyz"

    main(input_file)