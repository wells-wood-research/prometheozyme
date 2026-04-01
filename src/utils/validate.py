import os
from wsgiref import validate
import numpy as np
import re
from utils import structure
import pandas as pd

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

def evaluate_distance(conformation_coords, atom1_idx, atom2_idx):
    coord1 = conformation_coords[atom1_idx]
    coord2 = conformation_coords[atom2_idx]
    distance = calculate_distance(coord1, coord2)
    return distance

def evaluate_angle(conformation_coords, atom1_idx, atom2_idx, atom3_idx):
    coord1 = conformation_coords[atom1_idx]
    coord2 = conformation_coords[atom2_idx]
    coord3 = conformation_coords[atom3_idx]
    
    vec1 = coord1 - coord2
    vec2 = coord3 - coord2
    
    cos_angle = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
    angle_rad = np.arccos(np.clip(cos_angle, -1.0, 1.0))  # Clip to avoid numerical issues
    angle = np.degrees(angle_rad)
    return angle

def parse_energy_comment(comment):
    eopt = einter = None
    if comment:
        eopt_match = re.search(r"Eopt=(-?\d+\.\d+)", comment)
        einter_match = re.search(r"Einter=(-?\d+\.\d+)", comment)
        if eopt_match: eopt = float(eopt_match.group(1))
        if einter_match: einter = float(einter_match.group(1))
    return eopt, einter

def evaluate_restraints(coords, restraints):
    allOk = True
    for restr in restraints:
        atoms = restr.connectionsOpt

        if restr.type == "distance":
            currentValue = evaluate_distance(coords, atoms[0], atoms[1])

        elif restr.type == "angle":
            currentValue = evaluate_angle(coords, atoms[0], atoms[1], atoms[2])

        else:
            raise ValueError(f"Unknown restraint property: {restr.type}")

        # decide keep vs scan
        isOk = abs(currentValue - restr.value) <= restr.tolerance
        if not isOk:
                allOk = False
                break
    return allOk

def extract_docker_results(multi_xyz_path, assembly_metadata_comment, logger=None):
    """
    - Split a multi-structure XYZ file.
    - Extract Eopt/Einter from the comment line.
    - Add assembly's metadata in the comment line.
    - Save each as a single XYZ file.
    - Return list of (Result, DataFrame).
    """
    base_dir = os.path.dirname(multi_xyz_path)
    base_name = os.path.basename(multi_xyz_path)

    # Reuse universal parser
    structures = structure.read_multi_xyz(multi_xyz_path, logger=logger)

    results = []

    for i, (atom_count, comment, coords, atom_types) in enumerate(structures):
        # TODO after optmisation there might be duplicated results - need to remove based on RMSD
        # use AMPAL? https://isambard-uob.github.io/ampal/ampal.html#ampal.base_ampal.BaseAmpal.rmsd

        # --- Parse energies using regex ---
        eopt = einter = None
        if comment:
            eopt, einter = parse_energy_comment(comment)
            new_comment = f"{assembly_metadata_comment}||Eopt={eopt}|Einter={einter}"

        # --- Create dataframe ---
        df = pd.DataFrame({
            "ELEMENT": atom_types,
            "X": coords[:, 0],
            "Y": coords[:, 1],
            "Z": coords[:, 2],
        })

        # --- Construct new output filename ---
        new_path = os.path.join(base_dir, f"dock.result.{i}.xyz")

        # --- Reuse write_xyz() ---
        structure.write_xyz(new_path, new_comment, coords, atom_types)

        results.append((new_path, eopt, einter, df))

    if logger:
        logger.info(f"Split {multi_xyz_path} into {len(results)} results.")

    if len(results) < 1:
        logger.warning("No docker result structures passed the restraint test for this host seed. This is not critical unless no subsequent results exist. To increase success rate, consider increasing tolerance or force on the relevant restraint.")
    
    return results

def read_score(filepath, logger=None):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    # Find the start of the results table
    start_index = None
    for i, line in enumerate(lines):
        if line.strip().startswith("mode |  affinity"):
            start_index = i + 3  # data starts 3 lines after the header line
            break

    # Parse the affinity column
    affinities = []
    if start_index is not None:
        for line in lines[start_index:]:
            parts = line.split()
            if len(parts) >= 2:
                try:
                    affinity = float(parts[1])  # Second column is "affinity"
                    affinities.append(affinity)
                except ValueError:
                    continue  # Skip lines that don't contain floats in expected place
    else:
        logger.error(f"No scores found in {filepath}!")

    # Convert to NumPy array or DataFrame
    affinity_array = np.array(affinities)

    return affinity_array

def append_scores(xyz_file, scores_file, logger=None):
    structures = structure.read_xyz(xyz_file, logger=logger)
    scores = read_score(scores_file, logger=logger)

    if len(structures) != len(scores):
        if logger:
            logger.error(f"Mismatch: {len(structures)} structures vs {len(scores)} scores in {xyz_file}")
        else:
            raise ValueError(f"Mismatch: {len(structures)} structures vs {len(scores)} scores")

    temp_output = xyz_file + ".tmp"

    with open(temp_output, 'w') as f:
        for i, (atom_count, comment, coordinates, atom_types) in enumerate(structures):
            score = scores[i]
            new_comment = f"{score:.7f}"
            structure.write_xyz(f, new_comment, coordinates, atom_types)

    # Replace original file only after successful write
    os.replace(temp_output, xyz_file)