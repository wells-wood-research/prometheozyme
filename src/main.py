from collections import defaultdict
import logging
import os
import subprocess
import shutil
import sys
from pathlib import Path
import rmsd
import copy

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

def build_comment(species):
    return ";".join(f"{k}:{v}" for k, v in species.items())


def select_relevant_flavours(host, guest, idx_to_flavour):
    return {
        idx_to_flavour[col]
        for motif in (host, guest)
        for col, _ in motif
    }


def build_flavour_map(host, guest, idx_to_flavour):
    combined = dict(host)
    combined.update(dict(guest))

    return {
        idx_to_flavour[col]: (site.molecule_id, site.atom_id)
        for col, site in combined.items()
    }


def select_relevant_restraints(config, flavours):
    return [
        r for r in config.restraints.values()
        if all(conn in flavours for conn in r.connections)
    ]

def build_docking_restraints(base_restraints, flavour_map, guest_site, config):
    docking = []

    for r in base_restraints:

        # ORCA docker only supports distance restraints
        if r.type != "distance":
            continue

        sites = [flavour_map[f] for f in r.connections]

        if len(sites) != 2:
            continue  # skip invalid

        (mol1, atom1), (mol2, atom2) = sites

        # Convert node IDs → XYZ indices
        idx1 = config.ingredients[mol1].atoms[atom1].atomId
        idx2 = config.ingredients[mol2].atoms[atom2].atomId

        # Identify host vs guest
        if mol1 == guest_site.molecule_id:
            guest_atom = idx1
            host_atom = idx2
        elif mol2 == guest_site.molecule_id:
            guest_atom = idx2
            host_atom = idx1
        else:
            # skip host-host or guest-guest (not valid for docking)
            continue

        docking.append({
            "host_atom": host_atom,
            "guest_atom": guest_atom,
            "value": r.value,
            "tolerance": r.tolerance
        })

    return docking

def build_opt_restraints(base_restraints, flavour_map, host_node, guest_site, config):
    opt = []

    for r in base_restraints:
        new_r = copy.deepcopy(r)

        abs_atoms = []

        for fid in r.connections:
            mol_id, atom_id = flavour_map[fid]

            # Convert node ID → XYZ index
            local_idx = config.ingredients[mol_id].atoms[atom_id].atomId

            # Convert to assembly index
            if mol_id == guest_site.molecule_id:
                abs_idx = local_idx + host_node.n_atoms
            else:
                abs_idx = local_idx + host_node.species.get(mol_id, 0)

            abs_atoms.append(abs_idx)

        new_r.connectionsOpt = abs_atoms

        opt.append(new_r)

    return opt

# ----------------------------
# MAIN LOOP
# ----------------------------

def main(configPath):

    outdir, verbosity, rmsd_threshold, orca_params = parse.get_parameters()
    setup_logging(outdir, verbosity)

    config = parse.load_config(configPath)
    steps, dag, roots, rows = workflow.build_execution_plan(config.recipes)

    n_cols = len(config.flavours)
    lookUp = lookup.IDLookup()
    flavour_ids = list(config.flavours.keys())
    idx_to_flavour = dict(enumerate(flavour_ids))

    encoder = signature.SignatureEncoder(lookUp, n_cols)
    resolver = resolve.StructureFileResolver(encoder, config, cwd, outdir, rows)

    # -------------------------
    # INITIAL HOSTS
    # -------------------------
    for root_idx, root in enumerate(roots):
        hosts_by_motif = defaultdict(list)

        initial_files = resolver.resolve(root)

        for i, hf in enumerate(initial_files):
            _, _, coords, _ = structure.read_xyz(hf)

            mol_ids = {site.molecule_id for _, site in steps[0][0]}
            mol_id = next(iter(mol_ids))

            species = {mol_id: 0}

            hosts_by_motif[root].append(
                types.StepNode(
                    id=f"host_{i:03d}",
                    path=Path(hf),
                    step=0,
                    species=species,
                    n_atoms=len(coords)
                )
            )

        # -------------------------
        # MAIN LOOP - Track failed hosts for dependency skipping
        # -------------------------
        failed_host_motifs = set()
        
        for step_idx, (host, guest) in enumerate(steps):

            # Skip steps where the parent host motif has no valid candidates
            if host in failed_host_motifs:
                logging.warning(
                    f"Skipping step {step_idx} (parent motif has no valid hosts): "
                    f"{workflow.motif_to_row(host, n_cols)} -> {workflow.motif_to_row(guest, n_cols)}"
                )
                failed_host_motifs.add(guest)
                continue

            next_hosts = []

            guest_site, _ = workflow.determine_guest(host, guest)
            guest_file = structure.get_molec_xyz_path(
                cwd, config.ingredients[guest_site.molecule_id].filepath
            )

            guest_signature = encoder.motif_to_signature(guest)

            assembly_output_dir = structure.prep_assembly_dir(outdir / f"root_{root_idx:03d}", guest_signature)

            for host_idx, host_node in enumerate(hosts_by_motif[host]):

                host_dir = assembly_output_dir / f"host_{host_idx:03d}"
                dock_dir = host_dir / "dock"
                opt_dir = host_dir / "opt"
                valid_dir = host_dir / "valid"

                for d in [dock_dir, opt_dir, valid_dir]:
                    d.mkdir(parents=True, exist_ok=True)

                # -------------------------
                # PREP RESTRAINTS
                # -------------------------
                relevant_flavours = select_relevant_flavours(host, guest, idx_to_flavour)
                flavour_map = build_flavour_map(host, guest, idx_to_flavour)

                base_restraints = select_relevant_restraints(config, relevant_flavours)

                docking_restraints = build_docking_restraints(
                    base_restraints, flavour_map, guest_site, config
                )

                opt_restraints = build_opt_restraints(
                    base_restraints, flavour_map, host_node, guest_site, config
                )

                # -------------------------
                # PREP CHARGES AND MULTIPLICITIES
                # -------------------------
                host_charge = orca.calculate_charge([
                    config.ingredients[m].charge
                    for m in host_node.species.keys()
                ])

                host_multiplicity = orca.calculate_multiplicity([
                    config.ingredients[m].multiplicity
                    for m in host_node.species.keys()
                ])

                guest_charge = config.ingredients[guest_site.molecule_id].charge
                guest_multiplicity = config.ingredients[guest_site.molecule_id].multiplicity

                # -------------------------
                # DOCK
                # -------------------------
                dock_inp = orca.write_docking_input(
                    dock_dir,
                    host_idx,
                    host_node.path,
                    host_charge,
                    host_multiplicity,
                    guest_file,
                    guest_charge,
                    guest_multiplicity,
                    docking_restraints,
                    orca_params["qmMethod_dock"],
                    orca_params["strategy"],
                    orca_params["optLevel"],
                    orca_params["nOpt"],
                    orca_params["fixHost"],
                    orca_params["gridExtent"],
                    orca_params["nprocs"]
                )

                orca.run_orca(dock_inp, orca_params["orcapath"], timeout=None)

                dock_xyz = dock_inp.with_suffix(".docker.struc1.all.optimized.xyz")

                assembly_species = dict(host_node.species)
                assembly_species[guest_site.molecule_id] = host_node.n_atoms
                docked_results = validate.extract_docker_results(
                    dock_xyz,
                    build_comment(assembly_species),
                    logger=logging
                )

                # -------------------------
                # OPTIMISE
                # -------------------------
                for cand_idx, (path, _, _, df) in enumerate(docked_results):

                    candidate_dir = opt_dir / f"candidate_{cand_idx:03d}"
                    candidate_dir.mkdir(exist_ok=True)

                    coords = df[["X", "Y", "Z"]].to_numpy()

                    opt_inputs = orca.write_opt_input(
                        candidate_dir,
                        path,
                        coords,
                        orca.calculate_charge([host_charge, guest_charge]),
                        orca.calculate_multiplicity([host_multiplicity, guest_multiplicity]),
                        opt_restraints,
                        orca_params["qmMethod_opt"],
                        orca_params["nprocs"]
                    )

                    final_xyz = None

                    for inp in opt_inputs:
                        result = orca.run_orca(inp, orca_params["orcapath"], timeout=None)

                        if result != 0:
                            final_xyz = None
                            break  # stop chain immediately

                        final_xyz = inp.with_suffix(".xyz")

                    if final_xyz is None:
                        continue

                    # -------------------------
                    # VALIDATE
                    # -------------------------
                    coords = structure.read_xyz(final_xyz)[2]

                    if validate.evaluate_restraints(coords, opt_restraints):

                        assembly_comment = build_comment(assembly_species)

                        out_path = valid_dir / f"host_{cand_idx:03d}.xyz"
                        shutil.copy(final_xyz, out_path)
                        structure.write_xyz_comment(out_path, assembly_comment)

                        next_hosts.append(
                            types.StepNode(
                                id=f"{host_node.id}_c{cand_idx:03d}",
                                path=out_path,
                                step=step_idx + 1,
                                species=assembly_species,
                                n_atoms=len(coords)
                            )
                        )

            hosts_by_motif[guest].extend(next_hosts)

            if not next_hosts:
                logging.warning(
                    f"No valid hosts produced in step {step_idx}: "
                    f"{workflow.motif_to_row(host, n_cols)} -> {workflow.motif_to_row(guest, n_cols)}. "
                    f"Marking guest motif as failed and skipping dependent steps."
                )
                failed_host_motifs.add(guest)
                continue

            logging.info(
                f"Step {step_idx} complete: produced {len(hosts_by_motif[guest])} valid host(s)"
            )
    print("!" * 100)
    logging.info(f"Workflow complete. Processing results...")
    if failed_host_motifs:
        logging.warning(
            f"Failed to produce valid hosts for {len(failed_host_motifs)} motif(s). "
            f"Some recipe combinations could not be completed."
        )
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