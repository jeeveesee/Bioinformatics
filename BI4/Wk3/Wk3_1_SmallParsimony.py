##############################################################################################
# Molecular Evolution  - Bioinformatic IV Course from Coursera
#
# Week 3 - Small Parsimony (rooted tree)
# Code Challenge: Implement SmallParsimony to solve the Small Parsimony Problem.
# Input: An integer n followed by an adjacency list for a rooted binary tree with n leaves labeled by DNA strings.
# Output: The minimum parsimony score of this tree, followed by the adjacency list of a tree corresponding to
# labeling internal nodes by DNA strings in order to minimize the parsimony score of the tree.
# You may break ties however you like.
#
# Note: Remember to run SmallParsimony on each individual index of the strings at the leaves of the tree.
##############################################################################################
"""
Pseudocode (Compeau & Pevzner):

SmallParsimony(T, Character)
    for each node v in tree T
        Tag(v) ← 0
        if v is a leaf
            Tag(v) ← 1
            for each symbol k in the alphabet
                if Character(v) = k
                    sk(v) ← 0
                else
                    sk(v) ← ∞
    while there exist ripe nodes in T
        v ← a ripe node in T
        Tag(v) ← 1
        for each symbol k in the alphabet
            sk(v) ← minimumall symbols i {si(Daughter(v))+αi,k} + minimumall symbols j {sj(Son(v))+αj,k}
    return minimum over all symbols k {sk(v)}
"""

# Main code
#####################

from pathlib import Path as partho
from collections import defaultdict
from typing import Dict, List, Set, Tuple

def parse_adjacency_list(n: int, adjacency_list_str: str) -> Tuple[Dict[int, List[int]], Dict[int, int], Dict[int, str], int, Set[int]]:
    """
    Parse the adjacency list of a rooted binary tree with leaf DNA strings.

    Parameters:
        n (int): Number of leaves.
        adjacency_list_str (str): Multi-line string with parent->child mappings.

    Returns:
        adj (dict[int, list[int]]): Directed parent-to-children mapping.
        parents (dict[int, int]): Child-to-parent mapping.
        leaf_strings (dict[int, str]): Leaf node ID to DNA string mapping.
        root (int): Root node ID of the tree.
        all_nodes (set[int]): Set of all node IDs.
    """
    leaf_id = 0
    leaf_strings = {}
    adj = {}
    parents = {}
    all_nodes = set()
    
    for line in adjacency_list_str.strip().split("\n"):
        if "->" not in line:
            continue
        parent_str, child_str = line.strip().split("->")
        parent = int(parent_str)
        all_nodes.add(parent)
        
        if child_str.isdigit():
            child = int(child_str)
            all_nodes.add(child)
            if parent not in adj:
                adj[parent] = []
            adj[parent].append(child)
            parents[child] = parent
        else:
            child = leaf_id
            leaf_strings[child] = child_str
            leaf_id += 1
            all_nodes.add(child)
            if parent not in adj:
                adj[parent] = []
            adj[parent].append(child)
            parents[child] = parent
            
    # Find the root of the tree (node with no parent)
    root = None
    for node in all_nodes:
        if node not in parents:
            root = node
            break
            
    return adj, parents, leaf_strings, root, all_nodes


def get_post_order(root: int, adj: Dict[int, List[int]]) -> List[int]:
    """
    Compute a bottom-up post-order traversal of the tree.

    Parameters:
        root (int): Root node ID.
        adj (dict[int, list[int]]): Directed parent-to-children adjacency list.

    Returns:
        list[int]: Node IDs in post-order traversal order.
    """
    order = []
    def dfs(node):
        if node in adj:
            for child in adj[node]:
                dfs(child)
        order.append(node)
    dfs(root)
    return order

def small_parsimony(n: int, adjacency_list_str: str) -> Tuple[int, List[Tuple[str, str, int]]]:
    """
    Implement the Small Parsimony algorithm to reconstruct internal node sequences
    and find the minimum parsimony score.

    Parameters:
        n (int): Number of leaves.
        adjacency_list_str (str): Adjacency list representation of the tree.

    Returns:
        total_score (int): The minimum parsimony score of the tree.
        edges (list[tuple[str, str, int]]): List of directed edges with Hamming distance.
    """
    adj, parents, leaf_strings, root, all_nodes = parse_adjacency_list(n, adjacency_list_str)
    
    # Get sequence length from any leaf DNA string
    seq_len = len(next(iter(leaf_strings.values())))
    
    # Reconstructed characters for each node at each index
    reconstructed = {node: [None] * seq_len for node in all_nodes}
    for node in leaf_strings:
        reconstructed[node] = list(leaf_strings[node])
        
    # Standard alphabet order to match expected tie-breaking
    alphabet = 'ACTG'
    char_to_idx = {char: i for i, char in enumerate(alphabet)}
    idx_to_char = {i: char for i, char in enumerate(alphabet)}
    
    post_order = get_post_order(root, adj)
    total_score = 0
    
    # Run dynamic programming column by column
    for c in range(seq_len):
        s = {}
        for node in all_nodes:
            s[node] = [0] * 4
            
        # Initialize leaves
        for leaf_id, dna in leaf_strings.items():
            char = dna[c]
            target_idx = char_to_idx[char]
            for k in range(4):
                s[leaf_id][k] = 0 if k == target_idx else 999999
                
        # Bottom-up dynamic programming
        for node in post_order:
            if node in leaf_strings:
                continue
            left, right = adj[node]
            for k in range(4):
                min_left = min(s[left][i] + (0 if i == k else 1) for i in range(4))
                min_right = min(s[right][j] + (0 if j == k else 1) for j in range(4))
                s[node][k] = min_left + min_right
                
        # Top-down backtracking
        root_min_val = min(s[root])
        total_score += root_min_val
        
        # Pick root character (alphabetically first among minimizers)
        best_root_idx = s[root].index(root_min_val)
        reconstructed[root][c] = idx_to_char[best_root_idx]
        
        # Traverse top-down
        for node in reversed(post_order):
            if node == root or node in leaf_strings:
                continue
            parent = parents[node]
            parent_char = reconstructed[parent][c]
            parent_idx = char_to_idx[parent_char]
            
            # Find best character for current node among its minimizers
            node_min_val = min(s[node])
            s_minimizers = [k for k, val in enumerate(s[node]) if val == node_min_val]
            
            if parent_idx in s_minimizers:
                best_idx = parent_idx
            else:
                best_idx = s_minimizers[0]
                
            reconstructed[node][c] = idx_to_char[best_idx]
            
    # Format reconstructed strings
    reconstructed_strings = {node: "".join(chars) for node, chars in reconstructed.items()}
    
    # Reconstruct edges and distances
    edges = []
    for parent in adj:
        parent_str = reconstructed_strings[parent]
        for child in adj[parent]:
            child_str = reconstructed_strings[child]
            h_dist = sum(1 for a, b in zip(parent_str, child_str) if a != b)
            edges.append((parent_str, child_str, h_dist))
            edges.append((child_str, parent_str, h_dist))
            
    return total_score, edges


def formatterer(answer: Tuple[int, List[Tuple[str, str, int]]]) -> str:
    """
    Format the parsimony score and alphabetical list of edges for output.

    Parameters:
        answer (tuple[int, list[tuple[str, str, int]]]): Parsimony score and edge list.

    Returns:
        str: Multi-line formatted output string.
    """
    total_score, edges = answer
    lines = [str(total_score)]
    # Sort edges alphabetically to produce consistent, clean, and test-matching output
    sorted_edges = sorted(edges, key=lambda x: (x[0], x[1]))
    for u_str, v_str, h_dist in sorted_edges:
        lines.append(f"{u_str}->{v_str}:{h_dist}")
    return "\n".join(lines)


###########################################################################

if __name__ == "__main__":
# # Sample Input:
#     n = 4
#     adjacency_list_str = """
#     4->CAAATCCC
#     4->ATTGCGAC
#     5->CTGCGCTG
#     5->ATGGACGA
#     6->4
#     6->5"""

#     # Sample Output:
#     """16
#     ATTGCGAC->ATAGCCAC:2
#     ATAGACAA->ATAGCCAC:2
#     ATAGACAA->ATGGACTA:2
#     ATGGACGA->ATGGACTA:1
#     CTGCGCTG->ATGGACTA:4
#     ATGGACTA->CTGCGCTG:4
#     ATGGACTA->ATGGACGA:1
#     ATGGACTA->ATAGACAA:2
#     ATAGCCAC->CAAATCCC:5
#     ATAGCCAC->ATTGCGAC:2
#     ATAGCCAC->ATAGACAA:2
#     CAAATCCC->ATAGCCAC:5"""
#     answer = small_parsimony(n, adjacency_list_str)
#     print(formatterer(answer)) # Formatted answer

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

    answer = small_parsimony(n, adjacency_list_str)

    with open(current_dir / "Wk3_1_output.txt", "w") as output_file:
        output_file.write(formatterer(answer))

