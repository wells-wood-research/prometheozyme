import os
import numpy as np

def calculate_distance(coord1, coord2):
    """Calculate Euclidean distance between two 3D coordinates."""
    return np.sqrt(np.sum((coord1 - coord2) ** 2))

def calculate_center_of_mass(coordinates, indices, atom_types):
    """Calculate center of mass for given indices, weighted by atomic masses."""
    if not indices:
        return None
    
    # Atomic masses in atomic mass units (u) for common elements
    atomic_masses = {
        'H': 1.00794, 'He': 4.002602, 'C': 12.0107, 'N': 14.0067, 'O': 15.9994,
        'F': 18.998403, 'P': 30.973762, 'S': 32.065, 'Cl': 35.453,
        # Add more elements as needed
    }
    
    selected_coords = np.array([coordinates[i] for i in indices])
    masses = np.array([atomic_masses.get(atom_types[i], 1.0) for i in indices])  # Default to 1.0 if unknown
    
    # Weighted average: sum(mass * coord) / sum(mass)
    weighted_coords = selected_coords * masses[:, np.newaxis]
    center_of_mass = np.sum(weighted_coords, axis=0) / np.sum(masses)
    
    return center_of_mass

def evaluate_distance(conformation_coords, g_idx, h_idx, n_atoms_host):
    guest_coord = conformation_coords[g_idx + n_atoms_host]
    host_coord = conformation_coords[h_idx]
    distance = calculate_distance(guest_coord, host_coord)
    return distance