from itertools import combinations
from collections import defaultdict

def row_to_tuple(row):
    return tuple(entry.molecule_id for entry in row.values())

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
                    motif = dict(key)
                    motifs[frozenset(motif.items())] = set(row_ids)

    return motifs

def motif_cost(motif):
    # number of unique molecules (excluding "0" if desired)
    mols = {v for (_, v) in motif}
    return len(mols)


def motif_value(motif, covered_rows):
    support = len(covered_rows)
    specificity = len(motif)

    cost = motif_cost(motif)
    if cost == 0:
        return 0

    return support * specificity / cost

def is_submotif(a, b):
    # a ⊆ b
    return all(item in b for item in a)


def build_dag(motifs):
    dag = {m: set() for m in motifs}

    motif_list = list(motifs.keys())

    for i, a in enumerate(motif_list):
        for j, b in enumerate(motif_list):
            if i == j:
                continue

            if is_submotif(a, b):
                dag[a].add(b)

    return dag

def add_full_rows(rows, motifs):
    for i, row in enumerate(rows):
        motif = frozenset((idx, val) for idx, val in enumerate(row))
        motifs[motif] = {i}

# Select motifs that cover all rows with maximum reuse -
# weighted set cover problem, but small enough for branch & bound
def solve_optimal_plan(motifs):
    motif_items = list(motifs.items())

    best_plan = None
    best_score = float("-inf")

    def backtrack(i, covered, plan, score):
        nonlocal best_plan, best_score

        # all rows covered
        if len(covered) == TOTAL_ROWS:
            if score > best_score:
                best_score = score
                best_plan = plan[:]
            return

        if i >= len(motif_items):
            return

        motif, rows = motif_items[i]

        # Option 1: take motif
        new_rows = covered | rows
        gain = motif_value(motif, rows - covered)

        backtrack(
            i + 1,
            new_rows,
            plan + [motif],
            score + gain
        )

        # Option 2: skip
        backtrack(i + 1, covered, plan, score)

    TOTAL_ROWS = len(set().union(*motifs.values()))
    backtrack(0, set(), [], 0)

    return best_plan

def topo_sort(plan):
    plan_set = set(plan)

    edges = {m: set() for m in plan}

    for a in plan:
        for b in plan:
            if a != b and is_submotif(a, b):
                edges[a].add(b)

    visited = set()
    result = []

    def dfs(node):
        if node in visited:
            return
        visited.add(node)

        for nxt in edges[node]:
            dfs(nxt)

        result.append(node)

    for m in plan:
        dfs(m)

    return result[::-1]

def build_execution_plan(recipes):
    rows = [row_to_tuple(r) for r in recipes]

    motifs = generate_motifs(rows)
    add_full_rows(rows, motifs)

    dag = build_dag(motifs)

    optimal_plan = solve_optimal_plan(motifs)

    execution_order = topo_sort(optimal_plan)

    return execution_order