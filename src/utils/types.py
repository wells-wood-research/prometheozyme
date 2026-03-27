import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Literal

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
    filepath: str
    atoms: Dict[str, Atom]

# Inter-host/guest (new) restraints are introduced as BIAS in orca's %docker block 
# Other (old) restraints are written in orca's %geom block
# guest-guest restraints need to be considered in subsequent docker steps where guest has become host
Scope = Literal["host-guest", "host-host", "guest-guest"]

@dataclass
class Restraint:
    type: str
    value: float
    connections: List[str]
    connectionsTranslated: Optional[List[int]] = field(default_factory=list)
    currentValue: Optional[float] = None

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