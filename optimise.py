from drOrca import make_orca_input
import os

def calculate_charge(charges):
    return int(sum(charges))
    
def calculate_multiplicity(multiplicities):
    # Assume molecules are weakly or non-interacting
    spin = int(sum(m - 1 for m in multiplicities) / 2)
    multiplicity = 2*spin + 1
    return multiplicity

def optimise(arr, host_atom_count, logger):
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
    for guest_info in arr["guests_info"]:
        guest = guest_info["obj"]
        for constraint in guest.constraints:
            guestIdx, guestType, hostIdx, hostType, val = constraint
            # TODO hostType
            h_idx = hostIdx[0]
            if guestType == 'iter':
                for idx in guestIdx:
                    g_idx = host_atom_count + idx
                    atoms = [g_idx, h_idx]
                    geom["keep"].append({"atoms": atoms, "val": val})
            elif guestType == 'com':
                # TODO com
                g_idx = host_atom_count + guestIdx[0]
                atoms = [g_idx, hostIdx]
                geom["keep"].append({"atoms": atoms, "val": val})
            else:
                logger.error(f"guestType in constraints for role {arr.desc} cannot be {guestType}; only allowed options are 'iter' and 'com'!")

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