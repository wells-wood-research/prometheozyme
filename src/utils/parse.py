from typing import Dict, List
from utils.types import Atom, Config, Ingredient, Restraint, RecipeEntry, Recipes
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

def get_default_parameters():
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
    os.makedirs(outdir, exist_ok=True)

    return outdir, verbosity, rmsd_threshold, orca

def get_parameters(config):
        # Get miscellaneous parameters
    misc = config.get("misc", {})
    verbosity = misc.get("verbosity", ".")
    rmsd_threshold = misc.get("rmsd", 2.0)
    addTimestamp = misc.get("timestamp", True)
    project_name = misc.get("project_name", "test")
    workdir = misc.get("workdir", ".")

    # Get orca docker parameters
    orca = config.get("orca", {})

    # Setup output dir
    if addTimestamp:
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        project_name = f"{timestamp}_{project_name}"
    outdir = os.path.join(workdir, project_name)
    os.makedirs(outdir, exist_ok=True)

    # Get ingredients metadata
    ingredients = config.get("ingredients", [])

    return outdir, verbosity, rmsd_threshold, orca, ingredients

def get_cookbook(config):
    ingredients = config.get("ingredients", {})
    # dict[str, dict], keyed by ingredient id:
    # {
    #   id: {
    #     name: str,
    #     charge: int,
    #     multiplicity: int,
    #     file: str,
    #     atoms: dict[str, dict]  # keyed by atom id
    #       {
    #         atom_id: {
    #           elem: str,
    #           atomId: int
    #         }
    #       }
    #   }
    # }

    flavours = config.get("flavours", {})
    # dict[str, str]
    # { flavour_id: name }

    restraints = config.get("restraints", {})
    # dict[str, dict], keyed by restraint id:
    # {
    #   id: {
    #     type: str,
    #     value: float | None,
    #     connections: list[str]  # list of node ids
    #   }
    # }

    recipes = config.get("recipes", [])
    # list[dict[str, dict]]
    # each dict is keyed by flavour id:
    # {
    #   flavour_id: {
    #     molecule: { id: str },
    #     atom: { id: str }
    #   }
    # }

    return ingredients, flavours, restraints, recipes

def parse_ingredients(raw: dict) -> Dict[str, Ingredient]:
    ingredients: Dict[str, Ingredient] = {}

    for iid, data in raw.items():
        atoms_raw = data.get("atoms", {})

        atoms = {
            aid: Atom(
                elem=atom_data["elem"],
                atomId=atom_data["atomId"],
            )
            for aid, atom_data in atoms_raw.items()
        }

        ingredients[iid] = Ingredient(
            name=data["name"],
            charge=data["charge"],
            multiplicity=data["multiplicity"],
            filepath=data["file"],
            atoms=atoms,
        )

    return ingredients

def parse_flavours(raw: dict) -> Dict[str, str]:
    return dict(raw)  # already in desired format

def parse_restraints(raw: dict) -> Dict[str, Restraint]:
    restraints: Dict[str, Restraint] = {}

    for rid, data in raw.items():
        restraints[rid] = Restraint(
            type=data["type"],
            value=data.get("value"),
            connections=list(data.get("connections", [])),
        )

    return restraints

def parse_recipes(raw: list) -> Recipes:
    recipes: Recipes = []

    for recipe_dict in raw:
        parsed_recipe: Dict[str, RecipeEntry] = {}

        for flavour_id, entry in recipe_dict.items():
            parsed_recipe[flavour_id] = RecipeEntry(
                molecule_id=entry["molecule"]["id"],
                atom_id=entry["atom"]["id"],
            )

        recipes.append(parsed_recipe)

    return recipes

def parse_cookbook(ingredients_raw: dict, flavours_raw: dict, restraints_raw: dict, recipes_raw: list):   
    ingredients = parse_ingredients(ingredients_raw)
    flavours = parse_flavours(flavours_raw)
    restraints = parse_restraints(restraints_raw)
    recipes = parse_recipes(recipes_raw)

    return ingredients, flavours, restraints, recipes

def load_config(config_path: str) -> Config:
    raw_config = get_config(config_path)

    # parameters
    outdir, verbosity, rmsd_threshold, orca = get_default_parameters()
    # or: use get_parameters(raw_config) if switching later

    # cookbook
    ingredients_raw, flavours_raw, restraints_raw, recipes_raw = get_cookbook(raw_config)

    ingredients = parse_ingredients(ingredients_raw)
    flavours = parse_flavours(flavours_raw)
    restraints = parse_restraints(restraints_raw)
    recipes = parse_recipes(recipes_raw)

    return Config(
        ingredients=ingredients,
        flavours=flavours,
        restraints=restraints,
        recipes=recipes,
        outdir=outdir,
        verbosity=verbosity,
        rmsd_threshold=rmsd_threshold,
        orca=orca,
    )