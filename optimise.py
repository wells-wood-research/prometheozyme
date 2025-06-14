from drOrca import make_orca_input

def calculate_charge(charges):
    return int(sum(charges))
    
def calculate_multiplicity(multiplicities):
    # Assume molecules are weakly or non-interacting
    return int(sum(m - 1 for m in multiplicities) / 2)

def optimise(arr, path, host_atom_count):

    arrName = os.path.splitext(os.path.basename(arr["path"]))[0]
    orcaInputDir = os.path.join(os.path.dirname(arr["path"]), arrName) 
    os.makedirs(orcaInputDir, exist_ok=True)
    orcaInput = os.path.join(orcaInputDir, "opt.inp")

    title = f"""Initial optimisation of {arrName}:\n
                {arr["desc"]}"""
    
    qmMethod = "XTB2"
    inputFormat = "xyzfile"

    moleculeInfo = {"charge": calculate_charge([guest["obj"].charge for guest in arr["guests_info"]]),
                    "multiplicity": calculate_multiplicity([guest["obj"].multiplicity for guest in arr["guests_info"]])}

    parallelize = 8
    maxcore = 2500

    geom = {"keep": []}
    for guest_info in arr["guests_info"]:
        guest = guest_info["obj"]
        constraints = guest.constraints
        for constraint in guest.constraints:
            guestIdx, hostIdx, val = constraint
            guestIdx = host_atom_count + guestIdx
            atoms = [guestIdx, hostIdx]
            geom["keep"].append({"atoms": atoms, "val": val})

    make_orca_input(orcaInput = orcaInput,
                    title = title,
                    qmMethod = qmMethod,
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