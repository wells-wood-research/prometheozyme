import logging
from dataclasses import dataclass
from typing import Dict, List, Optional

########################
## LOGGING
########################

def get_logger(logger: Optional[logging.Logger] = None) -> logging.Logger:
    """Return a valid logger, creating one if not provided."""
    if logger is None:
        logger = logging.getLogger(__name__)
        if not logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')
            handler.setFormatter(formatter)
            logger.addHandler(handler)
        logger.setLevel(logging.INFO)
    return logger

########################
## CUSTOM CLASSES
########################

@dataclass
class Atom:
    elem: str
    atomId: int


@dataclass
class Ingredient:
    name: str
    charge: int
    multiplicity: int
    filepath: str  # ← add this
    atoms: Dict[str, Atom]


@dataclass
class Restraint:
    type: str
    value: Optional[float]
    connections: List[str]


@dataclass
class RecipeEntry:
    molecule_id: str
    atom_id: str


# recipes: list[dict[flavour_id, RecipeEntry]]
Recipes = List[Dict[str, RecipeEntry]]

@dataclass
class Config:
    ingredients: Dict[str, Ingredient]
    flavours: Dict[str, str]
    restraints: Dict[str, Restraint]
    recipes: Recipes

    outdir: str
    verbosity: str
    rmsd_threshold: float
    orca: dict

@dataclass(frozen=True)
class Site:
    molecule_id: str
    atom_id: str