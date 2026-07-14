##############################################################################################
# Molecular Evolution  - Bioinformatic IV Course from Coursera
#
# Week 2 - UPGMA (which stands for Unweighted Pair Group Method with Arithmetic Mean)
# Code Challenge: Implement UPGMA.
#
# Input: An integer n followed by a space separated n x n distance matrix.
# Output: An adjacency list for the ultrametric tree returned by UPGMA. Edge weights should be accurate to two decimal places
# (answers in the sample dataset below are provided to three decimal places).
# Note on formatting: The adjacency list must have consecutive integer node labels starting from 0.
# The n leaves must be labeled 0, 1, ..., n - 1 in order of their appearance in the distance matrix.
# Labels for internal nodes may be labeled in any order but must start from n and increase consecutively.
##############################################################################################
"""
Pseudocode (Compeau & Pevzner):

UPGMA(D, n)
    Clusters ← n single-element clusters labeled 1, ... , n
    construct a graph T with n isolated nodes labeled by single elements 1, ... , n
    for every node v in T
        Age(v) ← 0
    while there is more than one cluster
        find the two closest clusters Ci and Cj
        merge Ci and Cj into a new cluster Cnew with |Ci| + |Cj| elements
        add a new node labeled by cluster Cnew to T
        connect node Cnew to Ci and Cj by directed edges
        Age(Cnew) ← DCi, Cj / 2
        remove the rows and columns of D corresponding to Ci and Cj
        remove Ci and Cj from Clusters
        add a row/column to D for Cnew by computing D(Cnew, C) for each C in Clusters
        add Cnew to Clusters
    root ← the node in T corresponding to the remaining cluster
    for each edge (v, w) in T
        length of (v, w) ← Age(v) - Age(w)
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


def UPGMA(n, dist_matrix_str):
    """
    Build an ultrametric tree from a distance matrix using the UPGMA algorithm.

    Iteratively merges the two closest clusters, recording the new internal node
    and its age (half the merge distance). Edge weights are Age(parent) - Age(child).

    Parameters:
        n (int): number of leaf nodes
        dist_matrix_str (str): whitespace-separated n x n distance matrix string

    Returns:
        dict[int, list[tuple[int, float]]]: adjacency list mapping each node to
            a list of (neighbor, edge_weight) tuples (both directions stored)
    """
    raw = parse_dist_matrix(dist_matrix_str)

    # D[i][j] = distance between clusters i and j
    D = {i: {j: raw[i][j] for j in range(n)} for i in range(n)}

    cluster_size = {i: 1 for i in range(n)}
    age = {i: 0.0 for i in range(n)}
    adj = {i: [] for i in range(n)}

    clusters = list(range(n))
    next_label = n

    while len(clusters) > 1:
        # Find the closest pair of clusters
        min_dist = float('inf')
        ci, cj = -1, -1
        for a in range(len(clusters)):
            for b in range(a + 1, len(clusters)):
                u, v = clusters[a], clusters[b]
                if D[u][v] < min_dist:
                    min_dist = D[u][v]
                    ci, cj = u, v

        # Create new internal node
        cnew = next_label
        next_label += 1
        age[cnew] = min_dist / 2.0
        adj[cnew] = []

        # Connect cnew to ci and cj; edge weight = age difference
        w_ci = age[cnew] - age[ci]
        w_cj = age[cnew] - age[cj]
        adj[cnew].append((ci, w_ci))
        adj[cnew].append((cj, w_cj))
        adj[ci].append((cnew, w_ci))
        adj[cj].append((cnew, w_cj))

        # Compute distances from cnew to remaining clusters (weighted average)
        size_ci = cluster_size[ci]
        size_cj = cluster_size[cj]
        size_cnew = size_ci + size_cj
        cluster_size[cnew] = size_cnew
        D[cnew] = {}

        for c in clusters:
            if c == ci or c == cj:
                continue
            d_new = (size_ci * D[ci][c] + size_cj * D[cj][c]) / size_cnew
            D[cnew][c] = d_new
            D[c][cnew] = d_new

        D[cnew][cnew] = 0.0

        # Remove merged clusters from active set and distance table
        clusters.remove(ci)
        clusters.remove(cj)
        for c in clusters:
            D[c].pop(ci, None)
            D[c].pop(cj, None)
        del D[ci], D[cj]

        clusters.append(cnew)

    return adj


def formatterer(adj):
    """
    Format the UPGMA adjacency list into the required output string.

    Each directed edge is printed as "u->v:weight" with weight to three decimal places.
    Lines are sorted by source node, then by destination node.

    Parameters:
        adj (dict[int, list[tuple[int, float]]]): adjacency list from UPGMA

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
    # 0	20	17	11
    # 20	0	20	13
    # 17	20	0	10
    # 11	13	10	0"""

    # # Sample Output:
    # """0->5:7.000
    # 1->6:8.833
    # 2->4:5.000
    # 3->4:5.000
    # 4->2:5.000
    # 4->3:5.000
    # 4->5:2.000
    # 5->0:7.000
    # 5->4:2.000
    # 5->6:1.833
    # 6->5:1.833
    # 6->1:8.833"""
    # answer = UPGMA(n, dist_matrix_str)
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

    answer = UPGMA(n, dist_matrix_str)

    with open("Wk1_4_output.txt", "w") as output_file:
        output_file.write(formatterer(answer))
