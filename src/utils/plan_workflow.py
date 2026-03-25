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

def prune_motifs(motifs):
    motif_items = list(motifs.items())
    pruned = {}

    for i, (m1, r1) in enumerate(motif_items):
        dominated = False

        for j, (m2, r2) in enumerate(motif_items):
            if i == j:
                continue

            if set(m1).issubset(m2) and r1 == r2:
                dominated = True
                break

        if not dominated:
            pruned[m1] = r1

    return pruned
    
# Select motifs that cover all rows with maximum reuse -
# weighted set cover problem, but small enough for branch & bound
def plan_branch_and_bound(motifs):
    motif_items = list(motifs.items())

    # Precompute base values (max possible gain per motif)
    base_values = [
        motif_value(m, rows) for m, rows in motif_items
    ]

    # Sort motifs by descending value (important!)
    order = sorted(
        range(len(motif_items)),
        key=lambda i: base_values[i],
        reverse=True
    )

    motif_items = [motif_items[i] for i in order]
    base_values = [base_values[i] for i in order]

    # Precompute optimistic suffix sums for pruning
    suffix_max = [0] * (len(base_values) + 1)
    for i in range(len(base_values) - 1, -1, -1):
        suffix_max[i] = suffix_max[i + 1] + base_values[i]

    TOTAL_ROWS = len(set().union(*motifs.values()))

    best_score = float("-inf")
    best_plan = None

    # Memo: (index, frozenset(covered)) → best score seen
    memo = {}

    def backtrack(i, covered, plan, score):
        nonlocal best_score, best_plan

        covered_fs = frozenset(covered)

        # Memo pruning
        key = (i, covered_fs)
        if key in memo and memo[key] >= score:
            return
        memo[key] = score

        # All rows covered → valid solution
        if len(covered) == TOTAL_ROWS:
            if score > best_score:
                best_score = score
                best_plan = plan[:]
            return

        # No more motifs
        if i >= len(motif_items):
            return

        # Upper bound pruning
        optimistic = score + suffix_max[i]
        if optimistic <= best_score:
            return

        motif, rows = motif_items[i]

        # --- OPTION 1: TAKE ---
        new_rows = covered | rows
        gain = motif_value(motif, rows - covered)

        if gain > 0:  # skip useless additions
            backtrack(
                i + 1,
                new_rows,
                plan + [motif],
                score + gain
            )

        # --- OPTION 2: SKIP ---
        backtrack(i + 1, covered, plan, score)

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

    motifs = prune_motifs(motifs)   # 🔥 critical

    plan = plan_branch_and_bound(motifs)

    execution_order = topo_sort(plan)

    return execution_order