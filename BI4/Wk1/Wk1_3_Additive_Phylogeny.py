##############################################################################################
# Molecular Evolution  - Bioinformatic IV Course from Coursera
#
# Week 1 - Additive Phylogeny to solve the Distance-Based Phylogeny Problem
# Code Challenge: Implement AdditivePhylogeny to solve the Distance-Based Phylogeny Problem.
# Input: An integer n followed by a space-separated n x n distance matrix.
# Output: A weighted adjacency list for the simple tree fitting this matrix.
# Note on formatting: The adjacency list must have consecutive integer node labels starting from 0.
# The n leaves must be labeled 0, 1, ..., n - 1 in order of their appearance in the distance matrix.
# Labels for internal nodes may be labeled in any order but must start from n and increase consecutively.
##############################################################################################
"""
Pseudocode (Compeau & Pevzner):

AdditivePhylogeny(D)
    n ← number of rows in D
    if n = 2
        return the tree consisting of a single edge of length D1,2
    limbLength ← Limb(D, n)
    for j ← 1 to n - 1
        Dj,n ← Dj,n - limbLength
        Dn,j ← Dj,n
    (i, k) ← two leaves such that Di,k = Di,n + Dn,k
    x ← Di,n
    D ← D with row n and column n removed
    T ← AdditivePhylogeny(D)
    v ← the (potentially new) node in T at distance x from i on the path between i and k
    add leaf n back to T by creating a limb (v, n) of length limbLength
    return T
"""

import sys
from pathlib import Path as partho

sys.path.insert(0, str(partho(__file__).parent))
# parse_distance_matrix and limb_length reused from Wk1_2
# find_path has no equivalent in Wk1_1 or Wk1_2, so defined locally
from Wk1_2_Limb_Length import parse_distance_matrix, limb_length


# Main code


def find_path(tree, start, end):
    """
    Find the unique path between two nodes in a tree via DFS.

    Parameters:
        tree (dict[int, dict[int, int]]): Weighted adjacency dict.
        start (int): Source node.
        end (int): Destination node.

    Returns:
        list[int]: Ordered node labels from start to end.
    """
    stack = [(start, [start])]
    visited = set()
    while stack:
        node, path = stack.pop()
        if node == end:
            return path
        if node in visited:
            continue
        visited.add(node)
        for neighbor in tree.get(node, {}):
            if neighbor not in visited:
                stack.append((neighbor, path + [neighbor]))
    return None


def insert_node_on_path(tree, i, k, x, n_orig):
    """
    Return the node at distance x from i on the path i→k, creating it if necessary.

    If the walk lands on an existing node, return it. Otherwise split the edge at that
    point. New internal nodes are labelled max(tree, n_orig-1) + 1 so they never
    collide with leaf indices 0..n_orig-1.

    Parameters:
        tree (dict[int, dict[int, int]]): Weighted adjacency dict, modified in-place.
        i (int): Start node of the path.
        k (int): End node of the path.
        x (int): Distance from i at which to place the node.
        n_orig (int): Total number of leaves in the original problem.

    Returns:
        int: Node label at distance x from i on path i→k.
    """
    path = find_path(tree, i, k)
    dist = 0
    for idx in range(len(path) - 1):
        u, w = path[idx], path[idx + 1]
        edge_w = tree[u][w]
        if dist + edge_w == x:
            return w
        if dist + edge_w > x:
            v = max(max(tree), n_orig - 1) + 1
            rem = x - dist
            del tree[u][w]
            del tree[w][u]
            tree.setdefault(v, {})
            tree[u][v] = rem
            tree[v][u] = rem
            tree[v][w] = edge_w - rem
            tree[w][v] = edge_w - rem
            return v
        dist += edge_w
    return path[-1]


def additive_phylogeny(n, dist_matrix_str, _n_orig=None):
    """
    Reconstruct a simple weighted tree from an additive distance matrix.

    Parameters:
        n (int): Number of leaves in the current submatrix.
        dist_matrix_str (str): Multi-line tab-separated string of the n x n distance matrix.
        _n_orig (int): Total leaf count of the original problem; set automatically on first call.

    Returns:
        dict[int, dict[int, int]]: Weighted adjacency dict; leaves are 0..n-1,
        internal nodes start at n and increase consecutively.
    """
    if _n_orig is None:
        _n_orig = n
    d = parse_distance_matrix(dist_matrix_str)

    if n == 2:
        return {0: {1: d[0][1]}, 1: {0: d[0][1]}}

    limb = limb_length(n, n - 1, dist_matrix_str)

    for j in range(n - 1):
        d[j][n - 1] -= limb
        d[n - 1][j] = d[j][n - 1]

    # Find i, k such that D[i][k] = D[i][n-1] + D[n-1][k]
    i_star = k_star = -1
    for i in range(n - 1):
        for k in range(i + 1, n - 1):
            if d[i][k] == d[i][n - 1] + d[n - 1][k]:
                i_star, k_star = i, k
                break
        if i_star != -1:
            break

    x = d[i_star][n - 1]
    sub_d = [row[:-1] for row in d[:-1]]
    sub_str = '\n'.join('\t'.join(map(str, row)) for row in sub_d)
    T = additive_phylogeny(n - 1, sub_str, _n_orig)

    v = insert_node_on_path(T, i_star, k_star, x, _n_orig)
    T[v][n - 1] = limb
    T[n - 1] = {v: limb}
    return T


def formatterer(answer):
    """
    Format a weighted adjacency dict as a sorted 'u->v:w' adjacency list string.

    Parameters:
        answer (dict[int, dict[int, int]]): Weighted adjacency dict from additive_phylogeny.

    Returns:
        str: Each directed edge on its own line, sorted by source node then destination node.
    """
    lines = []
    for node in sorted(answer.keys()):
        for neighbor, weight in sorted(answer[node].items()):
            lines.append(f"{node}->{neighbor}:{weight}")
    return '\n'.join(lines)


###########################################################################

if __name__ == "__main__":
#     # Sample test
#     n = 4

#     dist_matrix_str = """0	13	21	22
# 13	0	12	13
# 21	12	0	13
# 22	13	13	0"""
#     # # Expected answer =
#     # 0->4:11
#     # 1->4:2
#     # 2->5:6
#     # 3->5:7
#     # 4->0:11
#     # 4->1:2
#     # 4->5:4
#     # 5->4:4
#     # 5->3:7
#     # 5->2:6
#     answer = additive_phylogeny(n, dist_matrix_str)
    # print(formatterer(answer)) # Formatted answer

    # From file

    # Get dataset
    from pathlib import Path as partho

    current_dir = partho(__file__).parent
    filename = input("Please enter the filename: ")
    file_path = current_dir / filename

    with open(file_path, "r") as file:
        data = file.read().strip()
        lines = data.split('\n')
        n = int(lines[0])
        dist_matrix = '\n'.join(lines[1:])

    answer = additive_phylogeny(n, dist_matrix)

    with open("Wk1_3_output.txt", "w") as output_file:
        output_file.write(formatterer(answer))
