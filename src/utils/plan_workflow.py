from typing import Dict, List
from collections import defaultdict
from utils.types import Recipes, RecipeEntry

def similarity(row_a: Dict[str, RecipeEntry], row_b: Dict[str, RecipeEntry]) -> int:
    return sum(1 for k in row_a if row_a[k] == row_b[k])


def unique_molecule_count(row: Dict[str, RecipeEntry]) -> int:
    return len({entry.molecule_id for entry in row.values()})


def order_group_by_similarity(group: List[Dict[str, RecipeEntry]]) -> List[Dict[str, RecipeEntry]]:
    """Greedy similarity ordering inside a group."""
    if not group:
        return []

    n = len(group)

    # Precompute similarity matrix
    sim = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            s = similarity(group[i], group[j])
            sim[i][j] = s
            sim[j][i] = s

    remaining = set(range(n))

    # Start with most "central" row
    start = max(range(n), key=lambda i: sum(sim[i]))
    order = [start]
    remaining.remove(start)

    while remaining:
        last = order[-1]
        next_idx = max(remaining, key=lambda j: sim[last][j])
        order.append(next_idx)
        remaining.remove(next_idx)

    return [group[i] for i in order]


def sort_recipes(recipes: Recipes) -> Recipes:
    # Step 1: group by unique molecule count
    groups = defaultdict(list)
    for row in recipes:
        groups[unique_molecule_count(row)].append(row)

    # Step 2: process groups in ascending order
    result = []
    for uniq_count in sorted(groups.keys()):
        group = groups[uniq_count]

        # Step 3: order within group by similarity
        ordered_group = order_group_by_similarity(group)

        result.extend(ordered_group)

    return result