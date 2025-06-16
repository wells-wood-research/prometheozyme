from drOrca import make_orca_input
import os
from utils import get_atom_count

def calculate_charge(charges):
    return int(sum(charges))
    
def calculate_multiplicity(multiplicities):
    # Assume molecules are weakly or non-interacting
    spin = int(sum(m - 1 for m in multiplicities) / 2)
    multiplicity = 2*spin + 1
    return multiplicity

def optimise(arr, host_atom_count, ingredient_map, logger):
    path = arr["path"]
    arrName = os.path.splitext(os.path.basename(path))[0]
    orcaInputDir = os.path.join(os.path.dirname(path), arrName) 
    os.makedirs(orcaInputDir, exist_ok=True)
    orcaInput = os.path.join(orcaInputDir, "opt.inp")

    title = f"Initial optimisation of {arrName}:\n # {arr['desc']}"
    
    qmMethod = "XTB2"
    method = "Opt"
    inputFormat = "xyzfile"

    moleculeInfo = {"charge": calculate_charge([guest["obj"].charge for guest in arr["guests_info"]]),
                    "multiplicity": calculate_multiplicity([guest["obj"].multiplicity for guest in arr["guests_info"]])}

    parallelize = 8
    maxcore = 2500

    # TODO com for guestType and all for hostType!
    geom = {"keep": []}
    accumulated_guest_atom_count = 0  # Initialize accumulated_guest_atom_count
    for guest_info in arr["guests_info"]:
        guest = guest_info["obj"]
        current_guest_atom_count = get_atom_count(ingredient_map[guest_info["obj"].name].path) # Get current guest's atom count
        current_guest_absolute_start_index = host_atom_count + accumulated_guest_atom_count # Calculate absolute start index for current guest
        for constraint in guest.constraints:
            guestIdx, guestType, hostIdx, hostType, val = constraint
            # TODO hostType
            h_idx = hostIdx[0]
            if guestType == 'iter':
                for idx in guestIdx:
                    g_idx_absolute = current_guest_absolute_start_index + idx # Calculate absolute g_idx
                    atoms = [g_idx_absolute, h_idx]
                    geom["keep"].append({"atoms": atoms, "val": val})
            elif guestType == 'com':
                # TODO com
                # Assuming guestIdx for 'com' is a single index relative to the current guest
                g_idx_absolute = current_guest_absolute_start_index + guestIdx[0] # Calculate absolute g_idx
                atoms = [g_idx_absolute, h_idx] # hostIdx might need adjustment if it's not absolute
                geom["keep"].append({"atoms": atoms, "val": val})
            else:
                logger.error(f"guestType in constraints for role {arr.desc} cannot be {guestType}; only allowed options are 'iter' and 'com'!")
        accumulated_guest_atom_count += current_guest_atom_count # Update accumulated_guest_atom_count

    make_orca_input(orcaInput = orcaInput,
                    title = title,
                    qmMethod = qmMethod,
                    method = method,
                    inputFormat = inputFormat,
                    inputFile = path,
                    moleculeInfo = moleculeInfo,
                    parallelize = parallelize,
                    maxcore = maxcore,
                    qmmm = None,
                    geom = geom,
                    neb = None,
                    scf = None,
                    docker = None)
    arr_optimised = orcaInput.replace(".inp", ".xyz")
    return arr_optimised

def pull_backbone(arr, host_atom_count, ingredient_map, logger, input_path):
    path = arr["path"]
    arrName = os.path.splitext(os.path.basename(path))[0]
    orcaInputDir = os.path.join(os.path.dirname(path), arrName) 
    os.makedirs(orcaInputDir, exist_ok=True)
    orcaInput = os.path.join(orcaInputDir, "pull.inp")

    title = f"Pulling backbones out of {arrName}:\n # {arr['desc']}"
    
    qmMethod = "XTB2"
    method = "Opt"
    inputFormat = "xyzfile"

    moleculeInfo = {"charge": calculate_charge([guest["obj"].charge for guest in arr["guests_info"]]),
                    "multiplicity": calculate_multiplicity([guest["obj"].multiplicity for guest in arr["guests_info"]])}

    parallelize = 8
    maxcore = 2500

    # TODO com for guestType and all for hostType!
    geom = {}
    geom["keep"] = []
    accumulated_guest_atom_count = 0  # Initialize accumulated_guest_atom_count
    for guest_info in arr["guests_info"]:
        guest = guest_info["obj"]
        current_guest_atom_count = get_atom_count(ingredient_map[guest_info["obj"].name].path) # Get current guest's atom count
        current_guest_absolute_start_index = host_atom_count + accumulated_guest_atom_count # Calculate absolute start index for current guest
        for constraint in guest.constraints:
            guestIdx, guestType, hostIdx, hostType, val = constraint
            # TODO hostType
            h_idx = hostIdx[0]
            g_idx_absolute = current_guest_absolute_start_index + guestIdx[0] # Calculate absolute g_idx
            atoms = [g_idx_absolute, h_idx] # hostIdx might need adjustment if it's not absolute
            geom["keep"].append({"atoms": atoms, "val": val})
        accumulated_guest_atom_count += current_guest_atom_count # Update accumulated_guest_atom_count
    geom["pots"] = []
    accumulated_guest_atom_count = 0  # Initialize accumulated_guest_atom_count
    for guest_info in arr["guests_info"]:
        guest = guest_info["obj"]
        current_guest_atom_count = get_atom_count(ingredient_map[guest_info["obj"].name].path) # Get current guest's atom count
        current_guest_absolute_start_index = host_atom_count + accumulated_guest_atom_count # Calculate absolute start index for current guest
        for constraint in guest.constraints:
            guestIdx, guestType, hostIdx, hostType, val = constraint
            val = 0.5 # Potentials have a constant force
            # TODO hostType
            h_idx = hostIdx[0]
            g_idx_absolute = current_guest_absolute_start_index # Calculate absolute g_idx
            atoms = [g_idx_absolute, h_idx] # hostIdx might need adjustment if it's not absolute
            geom["pots"].append({"atoms": atoms, "val": val})
        accumulated_guest_atom_count += current_guest_atom_count # Update accumulated_guest_atom_count

    make_orca_input(orcaInput = orcaInput,
                    title = title,
                    qmMethod = qmMethod,
                    method = method,
                    inputFormat = inputFormat,
                    inputFile = path,
                    moleculeInfo = moleculeInfo,
                    parallelize = parallelize,
                    maxcore = maxcore,
                    qmmm = None,
                    geom = geom,
                    neb = None,
                    scf = None,
                    docker = None)
    arr_optimised = orcaInput.replace(".inp", ".xyz")
    return arr_optimised