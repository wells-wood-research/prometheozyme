import yaml
import os
import sys
import datetime

def get_config(configPath):
    if not os.path.isfile(configPath):
        sys.exit(f"ERROR: File does not exist at {configPath}")
    with open(configPath, "r") as file:
        config = yaml.safe_load(file)
    return config

def get_standard_parameters():
    # misc = config.get("misc", {})
    verbosity = "debug" # misc.get("verbosity", ".")
    rmsd_threshold = 1.0 # misc.get("rmsd", 2.0)
    project_name = "Kemp_eliminase" # misc.get("project_name", "test")
    workdir = "/home/mchrnwsk/prometheozyme/runs" # misc.get("workdir", ".")

    # orca
    orca = {
        "orcapath": "/home/mchrnwsk/tools/orca_6_1_1_linux_x86-64_shared_openmpi418/orca",
        "qmMethod_dock": "XTB2",
        "qmMethod_opt": "XTB2",
        "strategy": "NORMAL",
        "optLevel": "sloppyopt",
        "nOpt": 5,
        "fixHost": True,
        "gridExtent": 15,
        "nprocs": 8
    }

    # Setup output dir
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    project_name = f"{timestamp}_{project_name}"
    outdir = os.path.join(workdir, project_name)

    return outdir, verbosity, rmsd_threshold, orca

def get_cookbook(config):
    ingredients = config.get("ingredients", []) # list of dicts, keys {id: string, name: string, atoms: list[dict{id: string, elem: string, atomId: int}]}
    flavours = config.get("flavours", []) # list of dicts, keys {id: string, name: string}
    restraints = config.get("restraints", []) # list of dicts, keys {id: string, type: string, value: float, connections: list[dict{id: string}]}
    recipes = config.get("recipes", []) # list of dicts, keys are flavour ids {id: dict{molecule: dict{id}}, atom: dict{id}}
    return ingredients, flavours, restraints, recipes

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Read configuration file.")
    parser.add_argument('--config', type=str, default=None, help="Path to the configuration YAML file")
    args = parser.parse_args()
    if not args.config:
        # TODO del static assignemnt when finished testing
        args.config = "/home/mchrnwsk/prometheozyme/graph-config.json"
    
    config = get_config(args.config) ## TODO it's a json not yaml - issue?
    outdir, verbosity, rmsd_threshold, orca = get_standard_parameters() ## TODO set parameters in config file and remove hardcoding
    ingredients, flavours, restraints, recipes = get_cookbook(config)
    print(ingredients)