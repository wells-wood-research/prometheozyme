import logging
import os
import subprocess
import shutil
import sys
from pathlib import Path
import rmsd
from collections import defaultdict
leftovers_dict = defaultdict(list)

from pathlib import Path
cwd = Path(__file__).resolve().parent

from utils import resolve, types, parse, workflow, structure, orca, validate, lookup, signature

def setup_logging(workdir, verbosity="info"):
    # Define a custom logging level called VERBOSE
    VERBOSE_LEVEL_NUM = 5
    logging.addLevelName(VERBOSE_LEVEL_NUM, "VERBOSE")
    def verbose_root(message, *args, **kwargs):
        logging.log(VERBOSE_LEVEL_NUM, message, *args, **kwargs)
    logging.verbose = verbose_root

    # Map string verbosity to logging levels
    level_map = {
        "verbose": VERBOSE_LEVEL_NUM,
        "debug": logging.DEBUG,
        "info": logging.INFO,
        "warning": logging.WARNING,
        "error": logging.ERROR,
        "critical": logging.CRITICAL,
    }
    log_level = level_map.get(verbosity.lower(), logging.INFO)
    
    # Remove existing handlers
    root_logger = logging.getLogger()
    for handler in list(root_logger.handlers):
        root_logger.removeHandler(handler)

    # Log to both console and file
    handlers = [logging.StreamHandler()]  # console
    handlers.append(logging.FileHandler(os.path.join(workdir, "log.txt")))  # file

    logging.basicConfig(
        format="%(asctime)s - %(levelname)s - %(message)s",
        level=log_level,
        handlers=handlers,
    )

# ----------------------------
# MAIN LOOP
# ----------------------------

def main(configPath):
    allOk = True
    
    # Read config file
    outdir, verbosity, rmsd_threshold, orca_params = parse.get_parameters() # TODO write parameters on export from frontend
    setup_logging(outdir, verbosity)
    logging.info(f"Cooking begins ({os.path.abspath(configPath)})")

    config = parse.load_config(configPath)
    steps, dag, roots, rows = workflow.build_execution_plan(config.recipes)
    workflow.print_execution_steps(steps, rows, logging)
    
    # Set up lookups etc
    n_cols = len(config.flavours)
    lookUp = lookup.IDLookup()
    flavour_ids = list(config.recipes[0].keys())
    idx_to_flavour = dict(enumerate(flavour_ids)) # flavour_to_idx = {v: k for k, v in idx_to_flavour.items()}
    def get_flav_ids(flavour_columns):
        return [idx_to_flavour[flav] for flav in flavour_columns]
    encoder = signature.SignatureEncoder(lookUp, n_cols)
    resolver = resolve.StructureFileResolver(encoder, config, cwd, outdir, rows)
    
    for host, guest in steps:
        # Find host file
        host_files = resolver.resolve(host)
        host_id, host_flav_cols = workflow.determine_host(host) # TODO must check that there is only one unique
            
        # Find guest file
        guest_id, guest_flav_cols = workflow.determine_guest(host, guest)
        guest_file = structure.get_molec_xyz_path(cwd, config.ingredients[guest_id.molecule_id].filepath)

        # Dock dir (for output of this assembly)
        guest_row = workflow.motif_to_row(guest, n_cols)
        guest_signature = encoder.motif_to_signature(guest)      
        assembly_output_dir = structure.prep_assembly_dir(outdir, guest_signature)

        # Used to map restrained connections to molecId, atomId sites
        flavour_map = lookup.build_flavour_map(guest_row, idx_to_flavour)
        # Find restraints to apply in this step
        host_flavours = get_flav_ids(host_flav_cols)
        guest_flavours = get_flav_ids(guest_flav_cols)
        flavours = host_flavours + guest_flavours
        relevant_restraints = [
            restr for restr in config.restraints.values()
            if all(connection in flavours for connection in restr.connections)
        ]

        host_assembly_n_atoms, host_assembly_metadata_comment, _, _ = structure.read_xyz(host_files[0]) # TODO host_files are multiple...
        host_assembly_species = {entry.split(':')[0]: int(entry.split(':')[1]) for entry in host_assembly_metadata_comment.split('||')[0].split(';')}
        # What do I need for the restraints?
        for restr in relevant_restraints:
            # what type: distance, angle, torsion
            sites = [flavour_map[fid] for fid in restr.connections]
            for site in sites:
                molecId = site[0]
                atomId = site[1]
                # whether they are between host-host, or host-guest, or guest-guest
                if molecId == guest_id.molecule_id:
                    scope =  "guest"
                    # or absolute assembly atom idxs to restrain simply add the number of atoms as guest is added to the end of the host file when making the assembly
                    atom_idx = config.ingredients[molecId].atoms[atomId].atomId
                else:
                    scope = "host"
                    # for absolute assembly atom idxs to restrain add firstAtomIdx of that molec in that assembly file described in xyz file's comment as {molecId: firstAtomIdx}
                    atom_idx = config.ingredients[molecId].atoms[atomId].atomId + host_assembly_species.get(molecId, 0)
                restr.connectionsDocking.append((scope, atom_idx))
            print(restr)
        host_charge = orca.calculate_charge([config.ingredients.get(molec, host_id.molecule_id).charge for molec in host_assembly_species.keys()])
        host_multiplicity = orca.calculate_multiplicity([config.ingredients.get(molec, host_id.molecule_id).multiplicity for molec in host_assembly_species.keys()])
        guest_charge = config.ingredients[guest_id.molecule_id].charge
        guest_multiplicity = config.ingredients[guest_id.molecule_id].multiplicity

        # assemble: dock, optimise, validate for success
        new_assembly_metadata_comment = f"{host_assembly_metadata_comment.split('||')[0]};{guest_id.molecule_id}:{host_assembly_n_atoms}"

        # dock - host frozen, distances between host and guest as bias
        docked_results = []
        for host_file in host_files:
            try:
                dock_inp_file_path = orca.write_docking_input(assembly_output_dir, host_file, host_charge, host_multiplicity, guest_file, guest_charge, guest_multiplicity, relevant_restraints, orca_params["qmMethod_dock"], orca_params["strategy"], orca_params["optLevel"], orca_params["nOpt"], orca_params["fixHost"], orca_params["gridExtent"], orca_params["nprocs"])
                logging.info(f"Docking input written at path: {dock_inp_file_path}")
                logging.info(f"Running docking...")
                result_dock = orca.run_orca(dock_inp_file_path, orca_params["orcapath"], timeout=None) # TODO orca needs to be a dataclass not dict
                logging.info(f"Docking complete. See details at path: {Path(dock_inp_file_path).with_suffix('.out')}")
                if result_dock != 0 :
                    logging.error(f"Docking returned status code {result_dock}. Skip to next theozyme...")
            except subprocess.CalledProcessError as e:
                logging.exception(f"Critical error during ORCA run: {e}")
            except Exception as e:
                logging.exception(f"Unexpected error in docking step: {e}")
            # process docking results - extract all possibilities from the multi-xyz output and assign assembly metadata comment
            docked_results.extend(validate.extract_docker_results(Path(dock_inp_file_path).with_suffix(".docker.struc1.all.optimized.xyz"), new_assembly_metadata_comment, logger=logging)) # TODO use optimised results as can control the number with preOpt
        
        # optimise
        # scan bonds, angles, dihedrals that have starting values from end of dock above tolerance away from target values in restraints - this should help with convergence and avoid issues with orca internal coordinates when angles approach linearity
        # keep other bonds, angles, dihedrals of host frozen
        # 3 scans at a time - batched chunks
        optimised_results = []
        candidate_id = 0
        for docked_result in docked_results: # each docked_result is (new_path, eopt, einter, df)
            try:
                file = Path(docked_result[0])
                coords = docked_result[3].loc[:, ["X", "Y", "Z"]].to_numpy()
                opt_output_dir = file.parent / f"candidate_{candidate_id}"
                opt_output_dir.mkdir()
                charge = orca.calculate_charge([host_charge, guest_charge])
                multiplicity = orca.calculate_multiplicity([host_multiplicity, guest_multiplicity])
                opt_restraints = [
                    types.Restraint(
                        type=r.type,
                        value=r.value,
                        connections=r.connections,
                        connectionsDocking=r.connectionsDocking,
                        connectionsOpt=[
                            b + host_assembly_n_atoms if a == "guest" else b
                            for a, b in r.connectionsDocking
                        ],
                        currentValue=r.currentValue
                    )
                    for r in relevant_restraints
                ]
                opt_inp_file_paths = orca.write_opt_input(opt_output_dir, file, coords, charge, multiplicity, opt_restraints, orca_params["qmMethod_opt"], orca_params["nprocs"])
                logging.info(f"Restrained optimisation input written at path(s): {opt_inp_file_paths}")
                logging.info(f"Running restrained optimisation...")
                for opt_inp_file_path in opt_inp_file_paths:
                    result_opt = orca.run_orca(opt_inp_file_path, orca_params["orcapath"], timeout=None)
                    logging.info(f"Optimisation complete: {opt_inp_file_path.with_suffix('.out')}")
                    if result_opt == 25: # issue with internal coordinates breaking when angle approaches linear
                        # check if opt.xyz exists - if yes, use, and maybe even proceed
                        logging.info("Optimisation failed due to issues with ORCA internal coordinates as the angle approaches linearity. Evaluating outputs for usability.")
                        if not os.path.exists(opt_inp_file_path.with_suffix(".xyz")):                    
                            # Rerun with COPT (in cartesian coordinates) - but then can't restrain anything
                            logging.info(f"Retrying optimisation in Cartesian coordinates")
                            copt_inp_file_path = orca.write_copt_input_path(opt_inp_file_path)
                            logging.info(f"Ooptimisation input rewritten in cartesian coordinates: {copt_inp_file_path}")
                            result_opt = orca.run_orca(copt_inp_file_path, orca_params["orcapath"], timeout=300) # TODO Put timeout in config
                            logging.info(f"Cartesian optimisation complete: {copt_inp_file_path.with_suffix('.out')}")
                    if result_opt != 0 :
                        logging.error(f"Optimisation returned status code {result_opt}. Skip to next theozyme...")
                        opt_inp_file_paths.remove(opt_inp_file_path)
            except Exception as e:
                logging.exception(f"Unexpected error in optimisation step: {e}")
            optimised_results.append((opt_inp_file_paths[-1].with_suffix(".xyz")))
            candidate_id += 1
        
        # check against restraints - passed become potential hosts for next steps
        for i, optimised_result in enumerate(optimised_results):
            coords = structure.read_xyz(optimised_result)[2]
            base_path = assembly_output_dir / "valid_result"
            extension = ".xyz"
            if validate.evaluate_restraints(coords, opt_restraints):
                dest_path = base_path.with_name(f"{base_path.stem}_{i}").with_suffix(extension)
                shutil.copy(optimised_result, dest_path)
                structure.write_xyz_comment(dest_path, new_assembly_metadata_comment)

    sys.exit(1) # TODO later add filtering based on rmsds
    # Create output directory
    results_dir = os.path.join(outdir, "results")
    specials_dir = os.path.join(outdir, "specials")
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(specials_dir, exist_ok=True)

    logging.info(f"Saving all successfully cooked theozymes to {results_dir}...")
    # Copy files with renamed names based on order
    for i, ing in enumerate(leftovers_list, start=1):
        new_name = f"result{i}"
        new_pathPDB = os.path.join(results_dir, f"{new_name}.pdb")
        # Copy and rename files
        shutil.copy(ing.pathPDB, new_pathPDB)
        # Log the results
        logging.info(
            f"path: {new_pathPDB}, eopt: {ing.eopt} (Eh), einter: {ing.einter} (kcal/mol)"
        )
    # Merge into one multi-frame PDB file for easier analysis
    result_paths = [os.path.join(results_dir, x) for x in os.listdir(results_dir) if x.lower().endswith(".pdb")]
    write_multi_pdb(result_paths, os.path.join(results_dir, "merged.pdb"))

    # Reduce on RMSD
    logging.info(f"Saving diverse theozymes with RMSD >= {rmsd_threshold} to {specials_dir}...")
    prod_count = 1
    dish_count = 1
    for meal in leftovers_dict.values():
        logging.debug(f"DEBUG processing dish {dish_count}")
        dish_count += 1
        diverse = []
        for ing in meal:
            logging.debug(f"DEBUG processing ingredient {ing_dish_signature(ing.df)}")
            # Compare only to already accepted structures
            isUnique = True
            min_rmsd = float('inf')
            for ref_ing in diverse:
                rmsd_val = float(rmsd.calculate_rmsd.main([ing.pathXYZ, ref_ing.pathXYZ]))
                if rmsd_val < min_rmsd:
                    min_rmsd = rmsd_val
                if rmsd_val < rmsd_threshold:
                    logging.debug(
                        f"Not saving {ing.name}: too similar to {ref_ing.name} (RMSD={rmsd_val:.3f})"
                    )
                    isUnique = False
                    break

            # If unique relative to ALL saved structures → keep it
            if isUnique:
                diverse.append(ing)
                new_name = f"result{prod_count}"
                new_pathPDB = os.path.join(specials_dir, f"{new_name}.pdb")
                shutil.copy(ing.pathPDB, new_pathPDB)
                logging.info(
                    f"Unique theozyme found: saved {ing.name} as {new_name}.pdb (Eopt={ing.eopt}, Einter={ing.einter}, min rmsd {min_rmsd})"
                )
                prod_count += 1

    # Merge into one multi-frame PDB file for easier analysis
    specials_paths = [os.path.join(specials_dir, x) for x in os.listdir(specials_dir) if x.lower().endswith(".pdb")]
    write_multi_pdb(specials_paths, os.path.join(specials_dir, "merged.pdb"))

    if allOk:
        logging.info("""
            Cooking complete - bon appetit!
            ⚝⭒٭⋆⚝⭒٭⋆⚝⭒٭⋆⚝⭒٭⋆⚝⭒٭⋆⚝⭒٭⋆
            """)
    else:
        logging.info("""
            Oh, crepe - something's gone wrong.
                     Try different settings?
            ☠︎︎⭒٭⋆☠︎︎⭒٭⋆☠︎︎⭒٭⋆☠︎︎⭒٭⋆☠︎︎⭒٭⋆☠︎︎⭒٭⋆
            """)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Process some ingredients and courses.")
    parser.add_argument('--graph_config', type=str, default=None, help="Path to the configuration YAML file generated by graph-based visual scripting website.")
    args = parser.parse_args()
    if not args.graph_config:
        # TODO del static assignemnt when finished testing
        args.graph_config = "/home/mchrnwsk/prometheozyme/graph-config.json"
    main(args.graph_config)