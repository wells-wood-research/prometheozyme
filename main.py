import yaml
import logging
from utils import Indices, Ingredient, Constraint, Role, update_guest_constraints, reduce_guests, print_reduced

# Configure logging to output to both console and file
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

# Console handler
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(formatter)
logger.addHandler(ch)

# File handler
fh = logging.FileHandler('main.log', mode='w')
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)
logger.addHandler(fh)

logger.info("Starting main script.")

configPath = "/home/mchrnwsk/theozymes/config.yaml"
# Load the YAML file
with open(configPath, "r") as file:
    config = yaml.safe_load(file)

ingredients = config.get("ingredients", [])
roles = config.get("roles", [])

logger.debug(f"""Ingredients are:
             {ingredients}""")
logger.debug(f"""Roles are:
             {roles}""")

ingredient_objects = []
ingredient_map = {}  # Map ingredient names to objects for role processing
for ing in ingredients:
    indices_dict = ing.get('indices', {})
    indices_obj = Indices(**indices_dict)  # Create Indices object with all YAML indices
    ingredient_obj = Ingredient(
        path=ing['path'],
        charge=ing['charge'],
        multiplicity=ing['multiplicity'],
        indices=indices_obj,
        name=ing['name']
    )
    ingredient_objects.append(ingredient_obj)
    ingredient_map[ing['name']] = ingredient_obj
logger.debug(f"""Ingredient map is:
            {ingredient_map}""")

# Find the substrate (host) ingredient, assumed to be 'sub'
host = next((ing for ing in ingredient_objects if ing.name == 'sub'), None)
if not host:
    print("Error: No 'sub' ingredient found for role host")
    exit(1)

# Create Role objects from YAML
role_objects = []
for role in roles:
    # Map candidates to Ingredient objects
    candidates = [ingredient_map[cand['name']] for cand in role['candidates']]
    
    # Handle constraints if present
    constraints_data = role.get('constraints', [])
    constraints = []
    for cons in constraints_data:
        constraint = Constraint(
            guestIdx=cons['guestIdx'],
            guestType=cons.get('guestType', 'iter'),
            hostIdx=cons['hostIdx'],
            hostType=cons.get('hostType', 'iter'),
            val=cons['val']
        )
        constraints.append(constraint)
    
    # Create Role object
    role_obj = Role(
        title=role['name'],
        priority=role['priority'],
        guests=candidates,
        host=host,
        constraints=constraints
    )
    role_objects.append(role_obj)
logger.debug(f"""Role objects are:
            {role_objects}""")

updated_guests = update_guest_constraints(role_objects)
logger.debug(f"""Updated guests are:
            {updated_guests}""")

unique_guests = reduce_guests(updated_guests)
logger.debug(f"""Unique guests are:
            {unique_guests}""")

print_reduced(updated_guests, unique_guests)

