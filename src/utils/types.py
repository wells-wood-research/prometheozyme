import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Literal
from pathlib import Path

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

# Allowed restraint types:
RestrType = Literal["distance", "angle", "dihedral"]

@dataclass
class Restraint:
    type: RestrType
    value: float
    tolerance: float = field(init=False)
    connections: List[str]
    connectionsDocking: Optional[List[int]] = field(default_factory=list)
    connectionsOpt: Optional[List[int]] = field(default_factory=list)
    currentValue: Optional[float] = None
    
    def __post_init__(self):
        self.tolerance = 0.5 if self.type == "distance" else 80.0 if self.type == "angle" else 100.0  # TODO assign in config

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

@dataclass
class StepNode:
    id: str
    path: Path
    step: int
    species: dict[str, int]   # {molecule_id: first_atom_index}
    n_atoms: int