##############################################################################################
# Molecular Evolution  - Bioinformatic IV Course from Coursera
#
# Week 3 - Nearest neighbors interchange heuristic (precursor to large parsimony)
# Code Challenge: Implement the nearest neighbor interchange heuristic for the Large Parsimony Problem.
#     Input: An integer n, followed by an adjacency list for an unrooted binary tree whose n leaves are labeled by DNA strings and whose internal nodes are labeled by integers.
#     Output: The parsimony score and unrooted labeled tree obtained after every step of the nearest neighbor interchange heuristic. Each step should be separated by a blank line.

# Note: Depending on how your code breaks ties, you may obtain a different solution than the one we provide.  As a result, the parsimony score at each step may vary.
##############################################################################################

# Pseudocode:
# NearestNeighborInterchange(Strings)
#      score ← ∞
#      generate an arbitrary unrooted binary tree Tree with |Strings| leaves
#      label the leaves of Tree by arbitrary strings from Strings
#      solve  the  Small Parsimony in an Unrooted Tree Problem for Tree
#      label the internal nodes of Tree according to a most parsimonious labeling
#      newScore ← the parsimony score of Tree
#      newTree ← Tree
#      while newScore < score
#           score ← newScore
#           Tree ← newTree
#           for each internal edge e in Tree
#                for each nearest neighbor NeighborTree of Tree with respect to the edge e
#                     solve the Small Parsimony in an Unrooted Tree Problem for NeighborTree
#                     neighborScore ← the minimum parsimony score of NeighborTree
#                     if neighborScore < newScore
#                          newScore ← neighborScore
#                          newTree ← NeighborTree
#      return newTree

# Main code here

import sys
from pathlib import Path as partho
from collections import defaultdict
from typing import Dict, List, Tuple

sys.path.insert(0, str(partho(__file__).parent))

from Wk3_2_SmallParsimony_Unrooted import small_parsimony_unrooted
from Wk3_3_NearestNeighbors_Tree import swap_neighbors
# Run file from BI3 folder as python -m Wk3_4_NearestNeighbors_Interchange (no py)


def parse_topology(adjacency_list_str: str) -> Dict[str, List[str]]:
    """
    Parse an unrooted binary tree's adjacency list into an undirected topology.

    Parameters:
        adjacency_list_str (str): Lines of the form "u->v" giving the tree's
            edges. Leaves are labeled by DNA strings; internal nodes are
            labeled by digit strings (e.g. "5", "6").

    Returns:
        dict[str, list[str]]: Undirected adjacency list keyed by node label,
            neighbors kept in first-seen order.
    """
    topology = defaultdict(list)
    for line in adjacency_list_str.strip().split("\n"):
        line = line.strip()
        if "->" not in line:
            continue
        u, v = line.split("->")
        if v not in topology[u]:
            topology[u].append(v)
    return dict(topology)


def get_internal_edges(topology: Dict[str, List[str]]) -> List[Tuple[str, str]]:
    """
    Find every internal edge of the tree, i.e. an edge joining two internal
    (digit-labeled) nodes rather than a leaf.

    Parameters:
        topology (dict[str, list[str]]): Undirected adjacency list of the tree.

    Returns:
        list[tuple[str, str]]: Each internal edge (a, b), listed once.
    """
    seen = set()
    edges = []
    for u in topology:
        if not u.isdigit():
            continue
        for v in topology[u]:
            if not v.isdigit():
                continue
            edge = frozenset((u, v))
            if edge not in seen:
                seen.add(edge)
                edges.append((u, v))
    return edges


def generate_nni_neighbors(topology: Dict[str, List[str]], a: str, b: str) -> List[Dict[str, List[str]]]:
    """
    Build the two nearest-neighbor rearrangements of the tree around edge (a, b).

    Parameters:
        topology (dict[str, list[str]]): Undirected adjacency list of the tree.
        a (str): One endpoint of the internal edge.
        b (str): Other endpoint of the internal edge.

    Returns:
        list[dict[str, list[str]]]: The two rearranged topologies.
    """
    w1, w2 = [node for node in topology[a] if node != b]
    x1, x2 = [node for node in topology[b] if node != a]

    return [
        swap_neighbors(topology, a, b, w2, x1),
        swap_neighbors(topology, a, b, w2, x2),
    ]


def score_topology(n: int, topology: Dict[str, List[str]]) -> Tuple[int, Dict[str, str], List[Tuple[str, str]]]:
    """
    Solve the Small Parsimony in an Unrooted Tree Problem for a given topology.

    Parameters:
        n (int): Number of leaves.
        topology (dict[str, list[str]]): Undirected adjacency list of the tree.

    Returns:
        tuple[int, dict[str, str], list[tuple[str, str]]]: Minimum parsimony
            score, reconstructed DNA strings keyed by topology label, and the
            directed edge list (both directions) for formatting.
    """
    adjacency_list_str = "\n".join(f"{u}->{v}" for u in topology for v in topology[u])
    return small_parsimony_unrooted(n, adjacency_list_str)


def nearest_neighbors_interchange(
    n: int, adjacency_list_str: str
) -> List[Tuple[int, Dict[str, str], List[Tuple[str, str]]]]:
    """
    Implement the Nearest Neighbor Interchange heuristic for the Large
    Parsimony Problem.

    Parameters:
        n (int): Number of leaves.
        adjacency_list_str (str): Adjacency list of a starting unrooted
            binary tree whose leaves are labeled by DNA strings and whose
            internal nodes are labeled by integers.

    Returns:
        list[tuple[int, dict[str, str], list[tuple[str, str]]]]: One entry
            per accepted step of the heuristic (the best neighbor tree found
            in each round that strictly improved on the round's starting
            score), each in the same shape returned by small_parsimony_unrooted
            (score, reconstructed strings, directed edges), ready for
            formatterer. The starting tree itself is never included, since it
            is not the result of an interchange.
    """
    topology = parse_topology(adjacency_list_str)

    score = float("inf")
    new_score, new_reconstructed, new_edges = score_topology(n, topology)
    new_tree = topology

    steps = []

    while new_score < score:
        score = new_score
        tree = new_tree

        for a, b in get_internal_edges(tree):
            for neighbor_tree in generate_nni_neighbors(tree, a, b):
                neighbor_score, neighbor_reconstructed, neighbor_edges = score_topology(n, neighbor_tree)
                if neighbor_score < new_score:
                    new_score = neighbor_score
                    new_tree = neighbor_tree
                    new_reconstructed = neighbor_reconstructed
                    new_edges = neighbor_edges

        if new_score < score:
            steps.append((new_score, new_reconstructed, new_edges))

    return steps


def formatterer(answer: List[Tuple[int, Dict[str, str], List[Tuple[str, str]]]]) -> str:
    """
    Format every accepted step of the heuristic into the expected output.

    Parameters:
        answer (list[tuple[int, dict[str, str], list[tuple[str, str]]]]): One
            entry per step, as returned by nearest_neighbors_interchange.

    Returns:
        str: The parsimony score and adjacency list for each step, with
            steps separated by a blank line, matching the expected output.
    """
    blocks = []
    for score, reconstructed_strings, original_edges in answer:
        lines = [str(score)]
        for u_str, v_str in original_edges:
            u_recon = reconstructed_strings[u_str]
            v_recon = reconstructed_strings[v_str]
            h_dist = sum(1 for a, b in zip(u_recon, v_recon) if a != b)
            lines.append(f"{u_recon}->{v_recon}:{h_dist}")
        blocks.append("\n".join(lines))

    return "\n\n".join(blocks)


###########################################################################

if __name__ == "__main__":
    # # Sample Input:
    # n = 5
    # adjacency_list_str = """
    # GCAGGGTA->5
    # TTTACGCG->5
    # CGACCTGA->6
    # GATTCCAC->6
    # 5->TTTACGCG
    # 5->GCAGGGTA
    # 5->7
    # TCCGTAGT->7
    # 7->5
    # 7->6
    # 7->TCCGTAGT
    # 6->GATTCCAC
    # 6->CGACCTGA
    # 6->7"""

    # # Sample Output:
    # """
    # 22
    # TCCGTAGT->TCAGCGGA:4
    # GATTCCAC->GAACCCGA:4
    # CGACCTGA->GAACCCGA:3
    # TTTACGCG->TCAGCGGA:5
    # TCAGCGGA->TTTACGCG:5
    # TCAGCGGA->GCAGCGGA:1
    # TCAGCGGA->TCCGTAGT:4
    # GCAGGGTA->GCAGCGGA:2
    # GCAGCGGA->GAACCCGA:3
    # GCAGCGGA->GCAGGGTA:2
    # GCAGCGGA->TCAGCGGA:1
    # GAACCCGA->GATTCCAC:4
    # GAACCCGA->CGACCTGA:3
    # GAACCCGA->GCAGCGGA:3

    # 21
    # TCCGTAGT->TCTGCGGA:4
    # GATTCCAC->GCTGCGGA:5
    # CGACCTGA->GCAGCGGA:4
    # TTTACGCG->TCTGCGGA:4
    # TCTGCGGA->TTTACGCG:4
    # TCTGCGGA->GCTGCGGA:1
    # TCTGCGGA->TCCGTAGT:4
    # GCAGGGTA->GCAGCGGA:2
    # GCTGCGGA->GCAGCGGA:1
    # GCTGCGGA->GATTCCAC:5
    # GCTGCGGA->TCTGCGGA:1
    # GCAGCGGA->CGACCTGA:4
    # GCAGCGGA->GCAGGGTA:2
    # GCAGCGGA->GCTGCGGA:1"""
    # answer = nearest_neighbors_interchange(n, adjacency_list_str)
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
        adjacency_list_str = '\n'.join(lines[1:])

    answer = nearest_neighbors_interchange(n, adjacency_list_str)

    with open(current_dir / "Wk3_4_output.txt", "w") as output_file:
        output_file.write(formatterer(answer))
