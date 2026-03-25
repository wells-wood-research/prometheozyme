from itertools import combinations
from collections import defaultdict, deque


# ----------------------------
# STEP 1 — Normalize rows
# ----------------------------
def row_to_tuple(row):
    return tuple(entry.molecule_id for entry in row.values())


# ----------------------------
# STEP 2 — Generate motifs
# ----------------------------
def generate_motifs(rows, min_support=2):
    n_cols = len(rows[0])
    motifs = {}

    for r in range(1, n_cols + 1):
        for cols in combinations(range(n_cols), r):
            groups = defaultdict(list)

            for i, row in enumerate(rows):
                key = tuple((c, row[c]) for c in cols)
                groups[key].append(i)

            for key, row_ids in groups.items():
                if len(row_ids) >= min_support:
                    motifs[frozenset(key)] = set(row_ids)

    return motifs


# ----------------------------
# STEP 3 — Add full rows
# ----------------------------
def add_full_rows(rows, motifs):
    for i, row in enumerate(rows):
        motif = frozenset((idx, val) for idx, val in enumerate(row))
        motifs[motif] = {i}


# ----------------------------
# STEP 4 — VALID edge rule
# One molecule per step (possibly multiple columns)
# ----------------------------
def is_valid_molecule_step(parent, child):
    parent_dict = dict(parent)
    child_dict = dict(child)

    # must be strict expansion
    if not parent_dict.items() < child_dict.items():
        return False

    # what was added
    added = {
        c: v for c, v in child_dict.items()
        if c not in parent_dict
    }

    if not added:
        return False

    # must involve only ONE molecule
    mols = set(added.values())
    if len(mols) != 1:
        return False

    mol = next(iter(mols))

    # find all columns where this molecule appears in the CHILD
    child_cols_for_mol = {
        c for c, v in child_dict.items() if v == mol
    }

    # ensure we added ALL of them at once
    return set(added.keys()) == child_cols_for_mol - set(parent_dict.keys())


# ----------------------------
# STEP 5 — Build DAG
# ----------------------------
def build_dag(motifs):
    dag = {m: set() for m in motifs}
    motif_list = list(motifs.keys())

    for a in motif_list:
        for b in motif_list:
            if a == b:
                continue

            if is_valid_molecule_step(a, b):
                dag[a].add(b)

    return dag


# ----------------------------
# STEP 6 — Find best root motif
# (max reuse, minimal size)
# ----------------------------
def find_root(motifs):
    best = None
    best_score = -1

    for m, rows in motifs.items():
        score = len(rows) / (len(m) + 1)

        if score > best_score:
            best_score = score
            best = m

    return best


# ----------------------------
# STEP 7 — BFS workflow
# ----------------------------
def bfs_workflow(dag, root):
    visited = set()
    queue = deque([root])

    workflow = []

    while queue:
        current = queue.popleft()

        if current in visited:
            continue

        visited.add(current)
        workflow.append(current)

        for child in dag[current]:
            queue.append(child)

    return workflow


# ----------------------------
# STEP 8 — Pretty printing
# ----------------------------
def motif_to_row(motif, n_cols):
    row = ["*"] * n_cols
    for c, v in motif:
        row[c] = v
    return row


def describe_step(parent, child):
    parent_dict = dict(parent)
    child_dict = dict(child)

    added = {
        c: v for c, v in child_dict.items()
        if c not in parent_dict
    }

    mol = next(iter(set(added.values())))
    cols = sorted(added.keys())

    return f"Add molecule {mol} at columns {cols}"


# ----------------------------
# STEP 9 — Build full plan
# ----------------------------
def build_execution_plan(recipes):
    rows = [row_to_tuple(r) for r in recipes]

    motifs = generate_motifs(rows)
    add_full_rows(rows, motifs)

    dag = build_dag(motifs)

    root = find_root(motifs)

    workflow = bfs_workflow(dag, root)

    return workflow, dag, root, rows


# ----------------------------
# STEP 10 — Print workflow nicely
# ----------------------------
def print_workflow(workflow, dag, root, rows):
    n_cols = len(rows[0])

    print("ROOT (start here):")
    print(motif_to_row(root, n_cols))
    print()

    step_num = 1

    for parent in workflow:
        for child in dag[parent]:
            if child in workflow:
                print(f"Step {step_num}:")
                print(" From:", motif_to_row(parent, n_cols))
                print(" To:  ", motif_to_row(child, n_cols))
                print(" Action:", describe_step(parent, child))
                print()
                step_num += 1