##############################################################################################
# Molecular Evolution  - Bioinformatic IV Course from Coursera
#
# Week 1 - Neighbor Joining Algorithm
# Code Challenge: Implement NeighborJoining.
#
# Input: An integer n, followed by an n x n distance matrix.
# Output: An adjacency list for the tree resulting from applying the neighbor-joining algorithm.
# Edge-weights should be accurate to two decimal places
# (they are provided to three decimal places in the sample output below).
#
# Note on formatting: The adjacency list must have consecutive integer node labels starting from 0.
# The n leaves must be labeled 0, 1, ..., n - 1 in order of their appearance in the distance matrix.
# Labels for internal nodes may be labeled in any order but must start from n and increase consecutively.
##############################################################################################
"""
Pseudocode (Compeau & Pevzner):

NeighborJoining(D)
    n ← number of rows in D
    if n = 2
        T ← tree consisting of a single edge of length D1,2
        return T
    D* ← neighbor-joining matrix constructed from the distance matrix D
    find elements i and j such that D*i,j is a minimum non-diagonal element of D*
    Δ ← (TotalDistanceD(i) - TotalDistanceD(j)) /(n - 2)
    limbLengthi ← (1/2)(Di,j + Δ)
    limbLengthj ← (1/2)(Di,j - Δ)
    add a new row/column m to D so that Dk,m = Dm,k = (1/2)(Dk,i + Dk,j - Di,j) for any k
    D ← D with rows i and j removed
    D ← D with columns i and j removed
    T ← NeighborJoining(D)
    add two new limbs (connecting node m with leaves i and j) to the tree T
    assign length limbLengthi to Limb(i)
    assign length limbLengthj to Limb(j)
    return T
"""

# Main code
#####################


def parse_dist_matrix(dist_matrix_str):
    """
    Parse a whitespace-separated distance matrix string into a 2D list of floats.
    Parameters:
        dist_matrix_str (str): multiline string of whitespace-separated values
    Returns:
        list[list[float]]: 2D distance matrix
    """
    lines = dist_matrix_str.strip().split('\n')
    return [list(map(float, line.split())) for line in lines]


def neighbor_joining(n, dist_matrix_str):
    """
    Build a phylogenetic tree from a distance matrix using the Neighbor Joining algorithm.
    Iteratively selects the pair (i, j) minimising the neighbor-joining criterion
    D*[i][j] = (n-2)*D[i][j] - TotalDist(i) - TotalDist(j), merges them into a
    new node m, and recurses on the reduced matrix. Limb-length connections for
    the removed pair are deferred and applied in reverse (post-order) so that the
    adjacency list insertion order mirrors the recursive formulation.
    Parameters:
        n (int): number of leaf nodes
        dist_matrix_str (str): whitespace-separated n x n distance matrix string
    Returns:
        dict[int, list[tuple[int, float]]]: adjacency list mapping each node to a
            list of (neighbor, edge_weight) tuples (both directions stored)
    """
    raw = parse_dist_matrix(dist_matrix_str)
    D = {i: {j: raw[i][j] for j in range(n)} for i in range(n)}
    adj = {i: [] for i in range(n)}

    clusters = list(range(n))
    next_label = n

    # Each entry: (m, ci, limb_i, cj, limb_j) deferred until after recursion
    deferred = []

    while len(clusters) > 2:
        m_size = len(clusters)
        totals = {i: sum(D[i][k] for k in clusters) for i in clusters}

        # Find (ci, cj) minimising D*[ci][cj]
        min_val = float('inf')
        ci, cj = -1, -1
        for a in range(len(clusters)):
            for b in range(a + 1, len(clusters)):
                u, v = clusters[a], clusters[b]
                d_star = (m_size - 2) * D[u][v] - totals[u] - totals[v]
                if d_star < min_val:
                    min_val = d_star
                    ci, cj = u, v

        delta = (totals[ci] - totals[cj]) / (m_size - 2)
        limb_i = (D[ci][cj] + delta) / 2
        limb_j = (D[ci][cj] - delta) / 2

        m = next_label
        next_label += 1
        adj[m] = []

        # Store limb connections to apply after inner recursion completes
        deferred.append((m, ci, limb_i, cj, limb_j))

        # Add new row/column for m: Dk,m = (Dk,i + Dk,j - Di,j) / 2
        D[m] = {}
        remaining = [c for c in clusters if c != ci and c != cj]
        for k in remaining:
            d_km = (D[k][ci] + D[k][cj] - D[ci][cj]) / 2
            D[m][k] = d_km
            D[k][m] = d_km
        D[m][m] = 0.0

        # Remove ci and cj from D and clusters
        clusters.remove(ci)
        clusters.remove(cj)
        for c in clusters:
            D[c].pop(ci, None)
            D[c].pop(cj, None)
        del D[ci], D[cj]

        clusters.append(m)

    # Base case: two clusters left — single edge
    u, v = clusters
    w = D[u][v]
    adj[u].append((v, w))
    adj[v].append((u, w))

    # Apply deferred limb connections in reverse (post-order) to mirror recursion
    for m, ci, limb_i, cj, limb_j in reversed(deferred):
        adj[m].append((ci, limb_i))
        adj[m].append((cj, limb_j))
        adj[ci].append((m, limb_i))
        adj[cj].append((m, limb_j))

    return adj


def formatterer(adj):
    """
    Format the neighbor-joining adjacency list into the required output string.
    Each directed edge is printed as "u->v:weight" with weight to three decimal
    places. Lines are sorted by source node, then by destination node.
    Parameters:
        adj (dict[int, list[tuple[int, float]]]): adjacency list from neighbor_joining
    Returns:
        str: formatted adjacency list string
    """
    lines = []
    for u in sorted(adj):
        for v, w in sorted(adj[u], key=lambda x: x[0]):
            lines.append(f"{u}->{v}:{w:.3f}")
    return '\n'.join(lines)


###########################################################################

if __name__ == "__main__":
# # Sample Input:
    # n = 4
    # dist_matrix_str = """
    # 0	23	27	20
    # 23	0	30	28
    # 27	30	0	30
    # 20	28	30	0"""

    # # Sample Output:
    # """0->4:8.000
    # 1->5:13.500
    # 2->5:16.500
    # 3->4:12.000
    # 4->5:2.000
    # 4->0:8.000
    # 4->3:12.000
    # 5->1:13.500
    # 5->2:16.500
    # 5->4:2.000"""
    # answer = neighbor_joining(n, dist_matrix_str)
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
        dist_matrix_str = '\n'.join(lines[1:])

    answer = neighbor_joining(n, dist_matrix_str)

    with open("Wk1_5_output.txt", "w") as output_file:
        output_file.write(formatterer(answer))

# # # Testa:
#     n = 4
#     dist_matrix_str = """
#     0     13	21	22
#     13    0	12	13
#     21	12	0	13
#     22	13	13	0"""

#     # # Sample Output:
# """0->4:11.000
# 1->4:2.000
# 2->5:6.000
# 3->5:7.000
# 4->0:11.000
# 4->1:2.000
# 4->5:4.000
# 5->2:6.000
# 5->3:7.000
# 5->4:4.000"""
#     answer = neighbor_joining(n, dist_matrix_str)
#     print(formatterer(answer)) # Formatted answer