##############################################################################################
# Molecular Evolution  - Bioinformatic IV Course from Coursera
#
# Week 1 - Distance matrices to evolutionary trees
# Distances Between Leaves Problem: Compute the distances between leaves in a weighted tree.
# Code Challenge: Solve the Distances Between Leaves Problem. The tree is given as an adjacency list of a graph whose leaves are integers between 0 and n - 1; the notation a->b:c means that node a is connected to node b by an edge of weight c. The matrix you return should be space-separated.
# Input:  An integer n followed by the adjacency list of a weighted tree with n leaves.
# Output: An n x n matrix (di,j), where di,j is the length of the path between leaves i and j.
##############################################################################################

# Pseudocode (Compeau & Pevzner):
#   DistancesBetweenLeaves(T, n):
#       for each pair of leaves i, j in {0..n-1}:
#           di,j ← length of the path connecting i and j in T
#       return (di,j)


def parse_adjacency_list(dist_matrix):
    """Parse a weighted adjacency list string into a nested dict graph.

    Parameters
    ----------
    dist_matrix : str
        Adjacency list where each line has the form 'a->b:c'.

    Returns
    -------
    dict[int, dict[int, int]]
        graph[u][v] = edge weight between u and v.
    """
    graph = {}
    for line in dist_matrix.strip().split('\n'):
        left, right = line.split('->')
        u = int(left)
        v_str, w_str = right.split(':')
        v, w = int(v_str), int(w_str)
        if u not in graph:
            graph[u] = {}
        graph[u][v] = w
    print(f"{graph=}")
    return graph



def distance_bw_leaves(n, dist_matrix):
    """Compute all pairwise distances between the n leaves of a weighted tree.

    Leaves are nodes labelled 0 through n-1. Internal nodes have higher labels.
    A DFS from each leaf accumulates the total edge weight along the unique
    tree path to every other reachable node.

    Parameters
    ----------
    n : int
        Number of leaves (labelled 0 … n-1).
    dist_matrix : str
        Adjacency list of the weighted tree (format 'a->b:c' per line).

    Returns
    -------
    list[list[int]]
        n x n matrix where entry [i][j] is the path length from leaf i to leaf j.
    """
    graph = parse_adjacency_list(dist_matrix)

    dist = [[0] * n for _ in range(n)]
    print(f"{dist=}")

    for leaf in range(n):
        # DFS — tree has no cycles, so first visit gives shortest (only) path
        visited = {}
        stack = [(leaf, 0)]
        print(f"{stack=}")
        while stack:
            node, d = stack.pop()
            print(f"{node=}")
            print(f"{d=}")
            if node in visited:
                continue
            visited[node] = d
            print(f"{visited=}")
            # dictionary.get(key, default=None)
            for neighbor, weight in graph.get(node, {}).items():
                print(f"{neighbor=}")
                print(f"{weight=}")
                if neighbor not in visited:
                    stack.append((neighbor, d + weight))
                    print(f"New stack, {stack=}")

        for other_leaf in range(n):
            dist[leaf][other_leaf] = visited.get(other_leaf, 0)

    return dist


def formatterer(answer):
    """Format an n x n distance matrix as tab-separated rows.

    Parameters
    ----------
    answer : list[list[int]]
        The distance matrix returned by distance_bw_leaves.

    Returns
    -------
    str
        Each row on its own line with values separated by tabs.
    """
    return '\n'.join('\t'.join(map(str, row)) for row in answer)


###########################################################################

if __name__ == "__main__":
    # Sample test
    n = 4
    dist_matrix = """0->4:11
1->4:2
2->5:6
3->5:7
4->0:11
4->1:2
4->5:4
5->4:4
5->3:7
5->2:6"""
    # Expected answer =
    # 0	13	21	22
    # 13  0	12	13
    # 21	12	0	13
    # 22	13	13	0
    answer = distance_bw_leaves(n, dist_matrix)
    print(formatterer(answer)) # Formatted answer

    # # From file

    # # Get dataset
    # from pathlib import Path as partho

    # current_dir = partho(__file__).parent
    # filename = input("Please enter the filename: ")
    # file_path = current_dir / filename

    # with open(file_path, "r") as file:
    #     data = file.read().strip()
    #     lines = data.split('\n')
    #     n = int(lines[0])
    #     dist_matrix = '\n'.join(lines[1:])

    # answer = distance_bw_leaves(n, dist_matrix)

    # with open("Wk1_1_output.txt", "w") as output_file:
    #     output_file.write(formatterer(answer))
