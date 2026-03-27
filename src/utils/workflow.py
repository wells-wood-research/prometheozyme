from itertools import combinations
from collections import defaultdict, deque
from utils.types import Site

# ----------------------------
# Utilities
# ----------------------------
def motif_to_row(motif, n_cols):
    row = [None] * n_cols
    for c, v in motif:
        row[c] = v
    return row

def describe_step(parent, child):
    parent_dict = dict(parent)
    child_dict = dict(child)
    added = {c: v for c, v in child_dict.items() if c not in parent_dict}
    if not added:
        return "No new molecules added"
    mol = next(iter(set(added.values())))
    cols = sorted(added.keys())
    return f"Add molecule {mol} at columns {cols}"

def determine_host(parent, child):
    parent_dict = dict(parent)
    added = {c: v for c, v in parent_dict.items()}
    if not added:
        return None, []
    mol = next(iter(set(added.values())))
    cols = sorted(added.keys())
    return mol, cols

def determine_guest(parent, child):
    parent_dict = dict(parent)
    child_dict = dict(child)
    added = {c: v for c, v in child_dict.items() if c not in parent_dict}
    if not added:
        return None, []
    mol = next(iter(set(added.values())))
    cols = sorted(added.keys())
    return mol, cols

def row_to_tuple(row):
    return tuple(
        Site(entry.molecule_id, entry.atom_id)
        for entry in row.values()
    )

def map_indices_to_column_ids(steps, recipes):
    """
    Map execution plan steps from numerical column indices to column IDs.

    steps: list of (parent, child) frozensets like {(0, '2'), (1, '5')}
    recipes: list of dicts (your config.recipes)
    """
    if not recipes:
        return steps

    # Assume all recipes have the same columns, take the first one
    col_ids = list(recipes[0].keys())  # preserves order
    index_to_colid = {i: cid for i, cid in enumerate(col_ids)}

    mapped_steps = []
    for parent, child in steps:
        parent_mapped = frozenset((index_to_colid[i], m) for i, m in parent)
        child_mapped  = frozenset((index_to_colid[i], m) for i, m in child)
        mapped_steps.append((parent_mapped, child_mapped))

    return mapped_steps

# ----------------------------
# Motif generation
# ----------------------------
def get_row_molecule_columns(row):
    mol_to_cols = {}
    for c, mol in enumerate(row):
        mol_to_cols.setdefault(mol, set()).add(c)
    return mol_to_cols

def is_complete_motif(motif, rows):
    motif_dict = dict(motif)
    for row in rows:
        mol_to_cols = get_row_molecule_columns(row)
        for mol in set(motif_dict.values()):
            cols_in_motif = {c for c, v in motif_dict.items() if v == mol}
            cols_in_row = mol_to_cols.get(mol, set())
            if cols_in_motif and cols_in_motif != cols_in_row:
                return False
    return True

def is_valid_motif(motif, mol_to_cols):
    motif_dict = dict(motif)
    for mol, cols in mol_to_cols.items():
        present = {c for c, v in motif_dict.items() if v == mol}
        if present and present != cols:
            return False
    return True

def generate_motifs(rows):
    n_cols = len(rows[0])
    motifs = {}
    for r in range(1, n_cols + 1):
        for cols in combinations(range(n_cols), r):
            groups = defaultdict(list)
            for i, row in enumerate(rows):
                key = tuple((c, row[c]) for c in cols)
                groups[key].append(i)
            for key, row_ids in groups.items():
                motif = frozenset(key)
                if not is_complete_motif(motif, [rows[i] for i in row_ids]):
                    continue
                motifs[motif] = set(row_ids)
    return motifs

def infer_molecule_groups(rows):
    mol_to_cols = {}
    for row in rows:
        for col, mol in enumerate(row):
            mol_to_cols.setdefault(mol, set()).add(col)
    return mol_to_cols

def add_full_rows(rows, motifs):
    for i, row in enumerate(rows):
        motif = frozenset((idx, val) for idx, val in enumerate(row))
        motifs[motif] = {i}

# ----------------------------
# DAG construction
# ----------------------------
def is_valid_molecule_step(parent, child, motifs, rows):
    parent_dict = dict(parent)
    child_dict = dict(child)
    if not parent_dict.items() < child_dict.items():
        return False
    added = {c: v for c, v in child_dict.items() if c not in parent_dict}
    if not added:
        return False
    added_mols = set(added.values())
    if len(added_mols) != 1:
        return False
    mol = next(iter(added_mols))
    row_ids = motifs[parent]
    allowed_cols = set()
    for rid in row_ids:
        row = rows[rid]
        for c, v in enumerate(row):
            if v == mol:
                allowed_cols.add(c)
    return set(added.keys()) == allowed_cols

def build_dag(motifs, rows):
    dag = {m: set() for m in motifs}
    motif_list = list(motifs.keys())
    for a in motif_list:
        for b in motif_list:
            if a != b and is_valid_molecule_step(a, b, motifs, rows):
                dag[a].add(b)
    return dag

def reduce_dag(dag):
    new_dag = {}
    for parent, children in dag.items():
        seen = set()
        new_children = set()
        for child in children:
            action = describe_step(parent, child)
            if action not in seen:
                seen.add(action)
                new_children.add(child)
        new_dag[parent] = new_children
    return new_dag

def add_missing_solution_layer(dag, motifs, rows):
    full_motifs = {
        frozenset((i, v) for i, v in enumerate(row))
        for row in rows
    }

    existing = set(motifs.keys())
    missing = full_motifs - existing

    for fm in missing:
        best_parent = max(existing, key=lambda m: len(m & fm))
        motifs[fm] = {i for i, row in enumerate(rows) if frozenset((c, v) for c, v in enumerate(row)) == fm}
        dag[best_parent].add(fm)
        dag[fm] = set()

# ----------------------------
# Root selection
# ----------------------------
def find_best_root(motifs, rows):
    full_row_ids = set(range(len(rows)))
    candidates = [m for m, rset in motifs.items() if rset == full_row_ids]
    if not candidates:
        return []
    def score(m):
        cols = len(m)
        mols = len(set(v for _, v in m))
        return (cols, -mols)
    best = max(candidates, key=score)
    return [best]

# ----------------------------
# BFS execution plan
# ----------------------------
def traverse_dag_unique(dag, roots):
    """
    Traverse DAG starting from roots and return a list of unique steps.
    Only record steps where new molecules are added.
    """
    seen = set(roots)
    queue = deque(roots)
    steps = []

    while queue:
        parent = queue.popleft()
        for child in dag.get(parent, []):
            if child in seen:
                continue
            # Only record if child adds at least one new molecule
            added = {c: v for c, v in child if c not in dict(parent)}
            if not added:
                seen.add(child)  # Still mark as seen to prevent future duplicates
                queue.append(child)
                continue

            seen.add(child)
            steps.append((parent, child))
            queue.append(child)

    return steps

# ----------------------------
# Visualization
# ----------------------------
def dag_levels(dag, roots):
    levels = defaultdict(list)
    visited = set()
    queue = deque([(r, 0) for r in roots])
    while queue:
        node, lvl = queue.popleft()
        if node in visited:
            continue
        visited.add(node)
        levels[lvl].append(node)
        for child in dag[node]:
            queue.append((child, lvl + 1))
    return levels

def visualize_dag(dag, roots, rows, logger=None):
    n_cols = len(rows[0])
    levels = dag_levels(dag, roots)
    logger.debug("\nDAG LEVELS:\n")
    for lvl in sorted(levels):
        logger.debug(f"Level {lvl}:")
        for node in levels[lvl]:
            logger.debug(f" {motif_to_row(node, n_cols)}")
        logger.debug("")
    logger.debug("\nEDGES:\n")
    for parent, children in dag.items():
        for child in children:
            logger.debug(f"{motif_to_row(parent, n_cols)} -> {motif_to_row(child, n_cols)} | {describe_step(parent, child)}")

# ----------------------------
# STEP 8 — Verification
# ----------------------------
def verify_all_reached(dag, roots, rows, logger=None):
    final_motifs = {frozenset((i, v) for i, v in enumerate(row)) for row in rows}
    visited = set()
    queue = deque(roots)
    while queue:
        node = queue.popleft()
        if node in visited:
            continue
        visited.add(node)
        for child in dag[node]:
            queue.append(child)
    missing = final_motifs - visited
    if missing:
        logger.warning("\nMissing solutions:")
        for m in missing:
            logger.warning(motif_to_row(m, len(rows[0])))

# ----------------------------
# Main pipeline
# ----------------------------
def build_execution_plan(recipes):
    rows = [row_to_tuple(r) for r in recipes]
    motifs = generate_motifs(rows)
    add_full_rows(rows, motifs)
    mol_to_cols = infer_molecule_groups(rows)
    motifs = {m: r for m, r in motifs.items() if is_valid_motif(m, mol_to_cols)}
    dag = build_dag(motifs, rows)
    dag = reduce_dag(dag)
    add_missing_solution_layer(dag, motifs, rows)
    roots = find_best_root(motifs, rows)
    steps = traverse_dag_unique(dag, roots)
    verify_all_reached(dag, roots, rows)
    return steps, dag, roots, rows

# ----------------------------
# Pretty printing
# ----------------------------
def print_execution_steps(steps, rows, logger=None):
    n_cols = len(rows[0])
    for i, (parent, child) in enumerate(steps, 1):
        logger.debug(f"""Step {i}: {describe_step(parent, child)}
                    From: {motif_to_row(parent, n_cols)}
                    To:   {motif_to_row(child, n_cols)}\n""")

def print_parallel_batches(dag, roots, logger=None):
    levels = dag_levels(dag, roots)
    logger.debug("\nPARALLEL EXECUTION PLAN:\n")
    for lvl in sorted(levels):
        if lvl == 0:
            continue
        logger.debug(f"Batch {lvl}:")
        for node in levels[lvl]:
            for parent, children in dag.items():
                if node in children:
                    logger.debug(f" {describe_step(parent, node)}")
        logger.debug("")