##############################################################################################
# Molecular Evolution  - Bioinformatic IV Course from Coursera
#
# Week 3 - Nearest neighbors of a tree problem (precursor to large parsimony)
# Code Challenge: Solve the Nearest Neighbors of a Tree Problem.
#   Input: Two internal nodes a and b specifying an edge e, followed by an adjacency list of an unrooted binary tree.
#   Output: Two adjacency lists representing the nearest neighbors of the tree with respect to e. Separate the adjacency lists with a blank line.
##############################################################################################

# Pseudocode:
#   1. Parse the edge (a, b) and the undirected adjacency list of the tree.
#   2. a has exactly two neighbors besides b: w1, w2.
#      b has exactly two neighbors besides a: x1, x2.
#   3. The two nearest-neighbor trees are formed by swapping w2 with x1, and w2 with x2,
#      i.e. detaching a subtree hanging off a and re-attaching it to b (and vice versa).
#   4. Format each resulting adjacency list, separating the two trees with a blank line.

# Main code here

from pathlib import Path as partho
from collections import defaultdict
from typing import Dict, List, Tuple


def parse_adjacency_list(adjacency_list_str: str) -> Tuple[str, str, Dict[str, List[str]]]:
    """
    Parse the edge (a, b) and adjacency list of an unrooted binary tree.

    Parameters:
        adjacency_list_str (str): First line "a b" naming the edge to rearrange
            around, followed by lines of the form "u->v" giving the tree's
            edges in both directions.

    Returns:
        a (str): Label of one endpoint of the edge to rearrange around.
        b (str): Label of the other endpoint of the edge to rearrange around.
        adjacency (dict[str, list[str]]): Undirected adjacency list keyed by
            node label, neighbors kept in first-seen order.
    """
    lines = [line.strip() for line in adjacency_list_str.strip().split("\n") if line.strip()]
    a, b = lines[0].split()

    adjacency = defaultdict(list)
    for line in lines[1:]:
        u, v = line.split("->")
        if v not in adjacency[u]:
            adjacency[u].append(v)

    return a, b, dict(adjacency)


def swap_neighbors(adjacency: Dict[str, List[str]], a: str, b: str, w: str, x: str) -> Dict[str, List[str]]:
    """
    Build one nearest-neighbor rearrangement of the tree around edge (a, b).

    Detaches w (currently a neighbor of a) and x (currently a neighbor of b),
    then reattaches w to b and x to a, producing one of the two nearest
    neighbors of the tree with respect to edge (a, b).

    Parameters:
        adjacency (dict[str, list[str]]): Undirected adjacency list of the original tree.
        a (str): One endpoint of the central edge.
        b (str): Other endpoint of the central edge.
        w (str): Neighbor of a to move over to b.
        x (str): Neighbor of b to move over to a.

    Returns:
        dict[str, list[str]]: Undirected adjacency list of the rearranged tree.
    """
    new_adjacency = {node: list(neighbors) for node, neighbors in adjacency.items()}

    new_adjacency[a][new_adjacency[a].index(w)] = x
    new_adjacency[b][new_adjacency[b].index(x)] = w
    new_adjacency[w][new_adjacency[w].index(a)] = b
    new_adjacency[x][new_adjacency[x].index(b)] = a

    return new_adjacency


def nearest_neighbors_tree(adjacency_list_str: str) -> List[Dict[str, List[str]]]:
    """
    Solve the Nearest Neighbors of a Tree Problem.

    Parameters:
        adjacency_list_str (str): First line "a b" naming the internal edge to
            rearrange around, followed by the adjacency list of an unrooted
            binary tree.

    Returns:
        list[dict[str, list[str]]]: The two nearest-neighbor trees with
            respect to edge (a, b), each as an undirected adjacency list.
    """
    a, b, adjacency = parse_adjacency_list(adjacency_list_str)

    w1, w2 = [node for node in adjacency[a] if node != b]
    x1, x2 = [node for node in adjacency[b] if node != a]

    tree_1 = swap_neighbors(adjacency, a, b, w2, x1)
    tree_2 = swap_neighbors(adjacency, a, b, w2, x2)

    return [tree_1, tree_2]


def formatterer(answer: List[Dict[str, List[str]]]) -> str:
    """
    Format the two nearest-neighbor trees as adjacency lists.

    Parameters:
        answer (list[dict[str, list[str]]]): The two nearest-neighbor trees,
            each as an undirected adjacency list.

    Returns:
        str: The two adjacency lists (edges sorted by node label), separated
            by a blank line, matching the expected output format.
    """
    tree_blocks = []
    for tree in answer:
        lines = []
        for node in sorted(tree, key=int):
            for neighbor in tree[node]:
                lines.append(f"{node}->{neighbor}")
        tree_blocks.append("\n".join(lines))

    return "\n\n".join(tree_blocks)


###########################################################################

if __name__ == "__main__":
    # # Sample Input:
    # adjacency_list_str = """
    # 5 4
    # 0->4
    # 4->0
    # 1->4
    # 4->1
    # 2->5
    # 5->2
    # 3->5
    # 5->3
    # 4->5
    # 5->4"""

    # # Sample Output:
    # """
    # 1->4
    # 0->5
    # 3->4
    # 2->5
    # 5->2
    # 5->4
    # 5->0
    # 4->1
    # 4->5
    # 4->3

    # 1->5
    # 0->4
    # 3->4
    # 2->5
    # 5->2
    # 5->4
    # 5->1
    # 4->0
    # 4->5
    # 4->3"""
    # answer = nearest_neighbors_tree(adjacency_list_str)
    # print(formatterer(answer)) # Formatted answer

    # From file

    # Get dataset
    from pathlib import Path as partho

    current_dir = partho(__file__).parent
    filename = input("Please enter the filename: ")
    file_path = current_dir / filename

    with open(file_path, "r") as file:
        adjacency_list_str = file.read().strip()

    answer = nearest_neighbors_tree(adjacency_list_str)

    with open(current_dir / "Wk3_3_output.txt", "w") as output_file:
        output_file.write(formatterer(answer))
