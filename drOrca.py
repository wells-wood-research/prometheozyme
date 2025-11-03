import re
from os import path as p
from subprocess import call
## CLEAN CODE ##
class FilePath:
    pass
class DirectoryPath:
    pass

# from https://github.com/wells-wood-research/qmmm/blob/559b8ddf88eaab55f138e4ad2bdc04c3016f1f2d/scripts/drOrca.py

def write_input(f, inputFormat, moleculeInfo, inputFile):
    qmCharge = moleculeInfo["charge"]
    qmMultiplicity = moleculeInfo["multiplicity"]

    f.write(f"*{inputFormat} {qmCharge} {qmMultiplicity} {inputFile}\n")
    f.write("\n")

def write_qmmm_block(f, qmmm):
    totalSystemInfo = qmmm.get("totalSystemInfo", {})
    if totalSystemInfo:
        totalSystemCharge = totalSystemInfo.get("charge", 0)
        totalSystemMultiplicity = totalSystemInfo.get("multiplicity", 1)

    parameterFile = qmmm.get("parameterFile", "")
    qmAtoms = qmmm.get("qmAtoms", "")
    UseInfoFromPDB = qmmm.get("UseInfoFromPDB", True)
    ExtendActiveRegion = qmmm.get("ExtendActiveRegion", "distance")
    DistExtendedActiveRegion = qmmm.get("DistExtendedActiveRegion", 1.0)
    
    f.write("%qmmm\n")
    f.write(f"{' '*4}Use_QM_InfoFromPDB {UseInfoFromPDB}\n")
    f.write(f"{' '*4}Use_Active_InfoFromPDB {UseInfoFromPDB}\n")
    if not UseInfoFromPDB and qmAtoms:
        f.write(f"{' '*4}QMAtoms {qmAtoms} end\n")
    if ExtendActiveRegion:
        f.write(f"{' '*4}ExtendActiveRegion {ExtendActiveRegion}\n")
        f.write(f"{' '*4}Dist_AtomsAroundOpt {DistExtendedActiveRegion}\n")
    f.write(f"{' '*4}Charge_Total {totalSystemCharge}\n")
    f.write(f"{' '*4}Mult_Total {totalSystemMultiplicity}\n")
    f.write(f"{' '*4}ORCAFFFilename \"{parameterFile}\"\n")
    f.write("end\n")
    f.write("\n")

def write_geometry_block(f, geom):
    geomBreak = geom.get("break", {})            
    geomKeep = geom.get("keep", {})
    geomScan = geom.get("scan", {})
    geomPots = geom.get("pots", {})

    f.write("%geom\n")
    if geomBreak:
        geomBreakAtoms = geomBreak.get("atoms", None)
        geomBreakAtoms = list(map(str, geomBreakAtoms))
        if len(geomBreakAtoms) == 2:
            geomBreakType = "B"
        f.write(f"{' '*4}modify_internal\n")
        f.write(f"{' '*8}{{ {geomBreakType} {' '.join(geomBreakAtoms)} R }}\n")
        f.write(f"{' '*4}end\n")
    if geomKeep:
        if isinstance(geomKeep, dict):
            geomKeepAtoms = geomKeep.get("atoms", None)
            geomKeepAtoms = list(map(str, geomKeepAtoms))
            if len(geomKeepAtoms) == 2:
                geomKeepType = "B"
            geomKeepVal = geomKeep.get("val", None)
            f.write(f"{' '*4}Constraints\n")
            f.write(f"{' '*8}{{ {geomKeepType} {' '.join(geomKeepAtoms)} {geomKeepVal} C }}\n")
            f.write(f"{' '*4}end\n")
        elif isinstance(geomKeep, list):
            f.write(f"{' '*4}Constraints\n")
            for constraint in geomKeep:
                geomKeepAtoms = constraint.get("atoms", None)
                geomKeepAtoms = list(map(str, geomKeepAtoms))
                if len(geomKeepAtoms) == 2:
                    geomKeepType = "B"
                geomKeepVal = constraint.get("val", None)
                f.write(f"{' '*8}{{ {geomKeepType} {' '.join(geomKeepAtoms)} {geomKeepVal} C }}\n")
            f.write(f"{' '*4}end\n")
    if geomScan:
        geomScanAtoms = geomScan.get("atoms", None)
        geomScanAtoms = list(map(str, geomScanAtoms))
        if len(geomScanAtoms) == 2:
            geomScanType = "B"
        scanValues = (geomScan.get("start", None), geomScan.get("end", None), geomScan.get("iter", None))
        f.write(f"{' '*4}Scan\n")
        f.write(f"{' '*8}{geomScanType} {' '.join(geomScanAtoms)} = {scanValues[0]}, {scanValues[1]}, {scanValues[2]}\n")
        f.write(f"{' '*4}end\n")
    if geomPots:
        if isinstance(geomPots, dict):
            geomPotsAtoms = geomPots.get("atoms", None)
            geomPotsAtoms = list(map(str, geomPotsAtoms))
            if len(geomPotsAtoms) == 2:
                geomPotsType = "B"
            geomPotsVal = geomPots.get("val", None)
            f.write(f"{' '*4}Potentials\n")
            f.write(f"{' '*8}{{ C {' '.join(geomPotsAtoms)} {geomPotsVal} }}\n")
            f.write(f"{' '*4}end\n")
        elif isinstance(geomPots, list):
            f.write(f"{' '*4}Potentials\n")
            for constraint in geomPots:
                geomPotsAtoms = constraint.get("atoms", None)
                geomPotsAtoms = list(map(str, geomPotsAtoms))
                if len(geomPotsAtoms) == 2:
                    geomPotsType = "B"
                geomPotsVal = constraint.get("val", None)
                f.write(f"{' '*8}{{ C {' '.join(geomPotsAtoms)} {geomPotsVal} }}\n")
            f.write(f"{' '*4}end\n")
    f.write("end\n")
    f.write("\n")

def write_neb_block(f, neb):
    nebMethod = neb.get("nebMethod", "NEB")
    nImages = neb.get("NImages", 7)
    printLevel = neb.get("printLevel", 1)
    abort = neb.get("abort", {})
    convergence = neb.get("convergence", {})
    springs = neb.get("springs", {})
    transRotDOF = neb.get("transRotDOF", {})
    initialPath = neb.get("initialPath", {})
    optimization = neb.get("optimization", {})
    idpp = initialPath.get("idpp", {})

    f.write("%neb\n")
    # Miscellaneous
    f.write(f"{' '*2}# Miscellaneous\n")
    f.write(f"{' '*2}NImages {nImages}\n")
    f.write(f"{' '*2}CI {nebMethod == 'NEB-CI'}\n")
    f.write(f"{' '*2}PrintLevel {printLevel}\n")
    f.write("\n")
    # Structures
    f.write(f"{' '*2}# Structures\n")
    f.write(f"{' '*2}Product_PDBFile \"{neb.get('productInput', '')}\"\n")
    f.write(f"{' '*2}PreOpt {neb.get('preOpt', False)}\n")
    f.write(f"{' '*2}Free_End {neb.get('freeEnd', False)}\n")
    f.write("\n")
    # Intermediates
    f.write(f"{' '*2}# Intermediates\n")
    f.write(f"{' '*2}NSteps_FoundIntermediate {abort.get('intermediateFoundSteps', 30)}\n")
    f.write("\n")
    # Convergence
    f.write(f"{' '*2}# Convergence\n")
    f.write(f"{' '*2}CheckSCFConv {abort.get('ifSCFunconverged', True)}\n")
    f.write(f"{' '*2}ConvType {convergence.get('convType', 'all')}\n")
    f.write(f"{' '*2}Tol_MaxFP_I {convergence.get('maxCompAtForcePerp2Path', 1.e-3)}\n")
    f.write(f"{' '*2}Tol_Scale {convergence.get('scaleTolRegVsCI', 10.0)}\n")
    f.write(f"{' '*2}Tol_MaxF_CI {convergence.get('maxCompAtForceOnCI', 2.e-3)}\n")
    f.write(f"{' '*2}Tol_RMSF_CI {convergence.get('rmsAtForceOnCI', 1.e-3)}\n")
    f.write(f"{' '*2}Tol_Turn_On_CI {convergence.get('minAtForceSwitchOnCI', 2.e-2)}\n")
    f.write("\n")
    # Springs
    f.write(f"{' '*2}# Springs\n")
    f.write(f"{' '*2}SpringType {springs.get('type', 'image')}\n")
    f.write(f"{' '*2}Energy_Weighted {springs.get('energyWeighted', True)}\n")
    f.write(f"{' '*2}SpringConst {springs.get('springConst', 0.01)}\n")
    f.write(f"{' '*2}SpringConst2 {springs.get('springConst2', 0.1)}\n")
    f.write(f"{' '*2}PerpSpring {springs.get('perpSpring', 'no')}\n")
    f.write(f"{' '*2}LLT_Cos {springs.get('lltCos', False)}\n")
    f.write("\n")
    # Translational and Rotational DOFs
    f.write(f"{' '*2}# Translational and rotational degrees of freedom\n")
    f.write(f"{' '*2}Quatern {transRotDOF.get('quatern', 'always')}\n")
    f.write(f"{' '*2}Fix_center {transRotDOF.get('fixCenter', True)}\n")
    f.write(f"{' '*2}Remove_extern_Force {transRotDOF.get('removeExternForce', True)}\n")
    f.write("\n")
    # Initial Path Generation
    f.write(f"{' '*2}# Initial path generation\n")
    f.write(f"{' '*2}Interpolation {initialPath.get('interpolation', 'IDPP')}\n")
    f.write(f"{' '*2}NPTS_Interpol {initialPath.get('nPtsInterpol', 10)}\n")
    f.write(f"{' '*2}Tangent {initialPath.get('tangent', 'improved')}\n")
    f.write(f"{' '*2}Prepare_Frags {initialPath.get('prepareFrags', True)}\n")
    f.write(f"{' '*2}Max_Frag_Dist {initialPath.get('maxFragDist', 3.5)}\n")
    f.write(f"{' '*2}Bond_Cutoff {initialPath.get('bondCutoff', 1.2)}\n")
    f.write(f"{' '*2}IDPP_NMax {idpp.get('nMax', 7000)}\n")
    f.write(f"{' '*2}IDPP_Tol_MaxF {idpp.get('tolMaxF', 0.01)}\n")
    f.write(f"{' '*2}IDPP_ksp {idpp.get('ksp', 1.0)}\n")
    f.write(f"{' '*2}IDPP_Alpha {idpp.get('alpha', 0.01)}\n")
    f.write(f"{' '*2}IDPP_MaxMove {idpp.get('maxMove', 0.05)}\n")
    f.write(f"{' '*2}IDPP_Debug {idpp.get('debug', True)}\n")
    f.write(f"{' '*2}IDPP_Quatern {idpp.get('quatern', True)}\n")
    f.write(f"{' '*2}IDPP_Dist_Interpolation {idpp.get('distInterpolation', 'Bilinear')}\n")
    f.write(f"{' '*2}IDPP_Bilinear_Partition {idpp.get('bilinearPartition', 0.5)}\n")
    f.write(f"{' '*2}SIDPP {idpp.get('sIdpp', False)}\n")
    f.write("\n")
    # Optimization Method
    f.write(f"{' '*2}# Optimisation method\n")
    f.write(f"{' '*2}Opt_Method {optimization.get('method', 'LBFGS')}\n")
    f.write(f"{' '*2}Maxmove {optimization.get('maxMove', 0.1)}\n")
    f.write(f"{' '*2}Stepsize {optimization.get('stepSize', 1.0)}\n")
    f.write(f"{' '*2}MaxIter {optimization.get('maxIter', 500)}\n")
    lbfgs = optimization.get("lbfgs", {})
    f.write(f"{' '*2}LBFGS_Memory {lbfgs.get('memory', 20)}\n")
    f.write(f"{' '*2}LBFGS_DR {lbfgs.get('dr', 1.e-3)}\n")
    f.write(f"{' '*2}LBFGS_Restart_On_Maxmove {lbfgs.get('restartOnMaxMove', True)}\n")
    f.write(f"{' '*2}LBFGS_Reparam_On_Restart {lbfgs.get('reparamOnRestart', False)}\n")
    f.write(f"{' '*2}LBFGS_Precondition {lbfgs.get('precondition', True)}\n")
    f.write("end\n")
    f.write("\n")

def write_scf_block(f, scf):
    f.write(f"%scf\n")
    f.write(f"{' '*4}KeepDens {scf.get('keepDens', False)}\n")
    f.write("end\n")
    f.write("\n")

def write_docker_block(f, docker):
    f.write(f"%DOCKER\n")
    f.write(f"{' '*4}GUEST \"{docker.get('guestPath', None)}\"\n")
    f.write(f"{' '*4}FIXHOST {docker.get('fixHost', True)}\n")
    if "bias" in docker.keys():
        biases = docker.get("bias", [])
        f.write(f"{' '*4}BIAS\n")
        for bias in biases:
            atoms = bias.get("atoms", [])
            val = bias.get("val", 0.0)
            force = bias.get("force", 100)
            f.write(f"{' '*8}{{ B {atoms[0]} {atoms[1]} {val} {force} }}\n")
        f.write(f"{' '*4}END\n")
    f.write("END\n")
    f.write("\n")

def write_simple_input(f, simpleInputLine, parallelize):
    for keyword in simpleInputLine:
        f.write(f"!{keyword}\n")
    if parallelize != 0:
        f.write(f"!PAL{parallelize}")
    f.write("\n")
    f.write("\n")

def write_title(f, title):
    x = len(title)
    f.write(f" # {'-'*(x+5)} #\n")
    f.write(f" # {title}{' '*5}\n")
    f.write(f" # {'-'*(x+5)} #\n")
    f.write("\n")

########################################

def modify_orca_input(input_file):
    # Read the input file
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    # Find the index of the '*' line (calculation start)
    star_index = next((i for i, line in enumerate(lines) if line.strip().startswith('*')), len(lines))
    
    # Search for '%geom' before star_index
    geom_start = next((i for i in range(star_index) if lines[i].strip().startswith('%geom')), None)
    
    if geom_start is not None:
        # Find the end of the %geom block
        open_blocks = 1
        i = geom_start + 1
        while i < star_index and open_blocks > 0:
            if lines[i].strip() == 'Constraints':
                open_blocks += 1
            elif lines[i].strip() == 'end':
                open_blocks -= 1
            i += 1
        geom_end = i - 1
        
        # Check if MaxIter is present, update or add it
        maxiter_found = False
        for j in range(geom_start + 1, geom_end):
            if lines[j].strip().startswith('MaxIter'):
                lines[j] = f'    MaxIter 500\n'
                maxiter_found = True
                break
        if not maxiter_found:
            lines.insert(geom_start + 1, f'    MaxIter 500\n')
            geom_end += 1
            if star_index < len(lines):
                star_index += 1
        
        # Find and modify the Constraints sub-block
        constraints_start = next((j for j in range(geom_start + 1, geom_end) if lines[j].strip() == 'Constraints'), None)
        if constraints_start is not None:
            # Find the end of Constraints
            open_subblocks = 1
            k = constraints_start + 1
            while k < geom_end and open_subblocks > 0:
                if lines[k].strip() == 'end':
                    open_subblocks -= 1
                k += 1
            constraints_end = k - 1
            
            # Modify constraint values
            constraint_pattern = r'\{\s*B\s+(\d+)\s+(\d+)\s+(\d+\.?\d*)\s*C\s*\}'
            for m in range(constraints_start + 1, constraints_end):
                match = re.search(constraint_pattern, lines[m])
                if match:
                    atom1, atom2, value = match.groups()
                    new_value = float(value) + 1.0
                    lines[m] = f'    {{ B {atom1} {atom2} {new_value:.1f} C }}\n'
    else:
        # No %geom block, add it before star_index
        new_geom = [
            f'%geom\n',
            f'    MaxIter 500\n',
            f'end\n'
        ]
        lines[star_index:star_index] = new_geom
    
    # Write back to the input file
    with open(input_file, 'w') as f:
        f.writelines(lines)
   
def make_orca_input(orcaInput: FilePath,
                    title: str,
                    simpleInputLine: str,
                    inputFormat: str,
                    inputFile: FilePath,
                    moleculeInfo: dict,
                    parallelize: int,
                    qmmm: dict,
                    geom: dict,
                    neb: dict,
                    scf: dict,
                    docker: dict
                    ) -> FilePath:
    
    if not simpleInputLine:
        simpleInputLine = "Opt" # Used to be "Opt", need to go back and fix what relies on it
    if neb:
        simpleInputLine.append(neb.get("nebMethod", "NEB"))
    if scf:
        if scf.get("keepDens", False):
            simpleInputLine.append("KeepDens")
    if qmmm:
        simpleInputLine.append("QMMM")

    if not orcaInput.endswith(".inp"):
        orcaInput += ".inp"
    with open(orcaInput, "w") as f:
        write_title(f, title)
        write_simple_input(f, simpleInputLine, parallelize)
        if geom:
            write_geometry_block(f, geom)
        if qmmm:
            write_qmmm_block(f, qmmm)
        if neb:
            write_neb_block(f, neb)
        if docker:
            write_docker_block(f, docker)
        write_input(f, inputFormat, moleculeInfo, inputFile)
    return orcaInput
