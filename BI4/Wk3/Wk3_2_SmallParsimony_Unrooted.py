##############################################################################################
# Molecular Evolution  - Bioinformatic IV Course from Coursera
#
# Week 3 - Small Parsimony in an unrooted tree problem
# Code Challenge: Solve the Small Parsimony in an Unrooted Tree Problem.
#    Input: An integer n followed by an adjacency list for an unrooted binary tree with n leaves labeled by DNA strings.
#    Output: The minimum parsimony score of this tree, followed by the adjacency list of the tree corresponding to labeling internal nodes 
#            by DNA strings in order to minimize the parsimony score of the tree.
##############################################################################################


from pathlib import Path as partho
from collections import defaultdict
from typing import Dict, List, Set, Tuple

def parse_unrooted_adjacency_list(adjacency_list_str: str) -> Tuple[defaultdict, Dict[str, int], Dict[int, str], List[Tuple[str, str]]]:
    """
    Parse the adjacency list of an unrooted binary tree.

    Parameters:
        adjacency_list_str (str): Multi-line string with undirected parent<->child mappings.

    Returns:
        undirected_adj (defaultdict[int, list[int]]): Undirected node ID adjacency list.
        node_to_id (dict[str, int]): Node label (string representation) to unique integer ID map.
        id_to_label (dict[int, str]): Unique integer ID to node label (string representation) map.
        original_edges (list[tuple[str, str]]): List of directed edges as they appeared in the input.
    """
    node_to_id = {}
    id_to_label = {}
    undirected_adj = defaultdict(list)
    original_edges = []
    
    for line in adjacency_list_str.strip().split("\n"):
        if "->" not in line:
            continue
        u_str, v_str = line.strip().split("->")
        original_edges.append((u_str, v_str))
        
        if u_str not in node_to_id:
            node_to_id[u_str] = len(node_to_id)
            id_to_label[node_to_id[u_str]] = u_str
        if v_str not in node_to_id:
            node_to_id[v_str] = len(node_to_id)
            id_to_label[node_to_id[v_str]] = v_str
            
        u = node_to_id[u_str]
        v = node_to_id[v_str]
        if v not in undirected_adj[u]:
            undirected_adj[u].append(v)
        if u not in undirected_adj[v]:
            undirected_adj[v].append(u)
            
    return undirected_adj, node_to_id, id_to_label, original_edges

def get_post_order(root: int, adj: Dict[int, List[int]]) -> List[int]:
    """
    Compute a bottom-up post-order traversal of the directed tree.

    Parameters:
        root (int): Root node ID.
        adj (dict[int, list[int]]): Directed parent-to-children adjacency list.

    Returns:
        list[int]: Node IDs in post-order traversal order.
    """
    order = []
    def dfs(node, parent):
        if node in adj:
            for child in adj[node]:
                if child != parent:
                    dfs(child, node)
        order.append(node)
    dfs(root, None)
    return order

def small_parsimony_unrooted(n: int, adjacency_list_str: str) -> Tuple[int, Dict[str, str], List[Tuple[str, str]]]:
    """
    Solve the Small Parsimony in an Unrooted Tree Problem.

    Parameters:
        n (int): Number of leaves.
        adjacency_list_str (str): Adjacency list representation of the unrooted tree.

    Returns:
        total_score (int): Minimum parsimony score.
        reconstructed_strings (dict[str, str]): Mapping of original node labels to their optimal DNA strings.
        original_edges (list[tuple[str, str]]): List of input edges for formatting in correct order.
    """
    undirected_adj, node_to_id, id_to_label, original_edges = parse_unrooted_adjacency_list(adjacency_list_str)
    
    # Identify leaf and internal nodes
    leaf_ids = [node_id for node_id in range(len(node_to_id)) if not id_to_label[node_id].isdigit()]
    
    leaf_strings = {leaf_id: id_to_label[leaf_id] for leaf_id in leaf_ids}
    seq_len = len(next(iter(leaf_strings.values())))
    
    # Select an arbitrary edge (u, v) to split and root the tree
    u, v = None, None
    for node in sorted(undirected_adj.keys()):
        if undirected_adj[node]:
            u = node
            v = undirected_adj[node][0]
            break
            
    # Remove edge (u, v) from undirected_adj
    undirected_adj[u].remove(v)
    undirected_adj[v].remove(u)
    
    # Root the tree at a new root node
    root_id = len(node_to_id)
    
    # Build directed tree
    directed_adj = defaultdict(list)
    parents = {}
    visited = {root_id}
    
    def build_directed_tree(curr, parent):
        visited.add(curr)
        parents[curr] = parent
        if parent is not None and parent != root_id:
            directed_adj[parent].append(curr)
        for neighbor in undirected_adj[curr]:
            if neighbor not in visited:
                build_directed_tree(neighbor, curr)
                
    # Root has two children: u and v
    directed_adj[root_id] = [u, v]
    parents[u] = root_id
    parents[v] = root_id
    
    build_directed_tree(u, root_id)
    build_directed_tree(v, root_id)
    
    all_nodes = set(range(root_id + 1))
    
    # Reconstructed characters for each node at each index
    reconstructed = {node: [None] * seq_len for node in all_nodes}
    for leaf_id in leaf_strings:
        reconstructed[leaf_id] = list(leaf_strings[leaf_id])
        
    alphabet = 'ACTG'
    char_to_idx = {char: i for i, char in enumerate(alphabet)}
    idx_to_char = {i: char for i, char in enumerate(alphabet)}
    
    post_order = get_post_order(root_id, directed_adj)
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
            left, right = directed_adj[node]
            for k in range(4):
                min_left = min(s[left][i] + (0 if i == k else 1) for i in range(4))
                min_right = min(s[right][j] + (0 if j == k else 1) for j in range(4))
                s[node][k] = min_left + min_right
                
        # Top-down backtracking
        root_min_val = min(s[root_id])
        total_score += root_min_val
        
        # Pick root character (alphabetically first among minimizers)
        best_root_idx = s[root_id].index(root_min_val)
        reconstructed[root_id][c] = idx_to_char[best_root_idx]
        
        # Traverse top-down
        for node in reversed(post_order):
            if node == root_id or node in leaf_strings:
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
            
    # Format reconstructed strings (excluding temporary root_id)
    reconstructed_strings = {id_to_label[node]: "".join(chars) for node, chars in reconstructed.items() if node != root_id}
    
    return total_score, reconstructed_strings, original_edges

def formatterer(answer: Tuple[int, Dict[str, str], List[Tuple[str, str]]]) -> str:
    """
    Format the parsimony score and list of reconstructed edges.

    Parameters:
        answer (tuple[int, dict[str, str], list[tuple[str, str]]]): Total score, reconstructed strings, and original edge pairs.

    Returns:
        str: Formatted output string matching sample output format exactly.
    """
    total_score, reconstructed_strings, original_edges = answer
    lines = [str(total_score)]
    
    for u_str, v_str in original_edges:
        u_recon = reconstructed_strings[u_str]
        v_recon = reconstructed_strings[v_str]
        h_dist = sum(1 for a, b in zip(u_recon, v_recon) if a != b)
        lines.append(f"{u_recon}->{v_recon}:{h_dist}")
        
    return "\n".join(lines)


###########################################################################

if __name__ == "__main__":
    # # Sample Input:
    # n = 4
    # adjacency_list_str = """
    # TCGGCCAA->4
    # 4->TCGGCCAA
    # CCTGGCTG->4
    # 4->CCTGGCTG
    # CACAGGAT->5
    # 5->CACAGGAT
    # TGAGTACC->5
    # 5->TGAGTACC
    # 4->5
    # 5->4"""

    # # Sample Output:
    # """17
    # TCGGCCAA->CCAGGCAC:4
    # CCTGGCTG->CCAGGCAC:3
    # TGAGTACC->CAAGGAAC:4
    # CCAGGCAC->CCTGGCTG:3
    # CCAGGCAC->CAAGGAAC:2
    # CCAGGCAC->TCGGCCAA:4
    # CACAGGAT->CAAGGAAC:4
    # CAAGGAAC->CACAGGAT:4
    # CAAGGAAC->TGAGTACC:4
    # CAAGGAAC->CCAGGCAC:2"""
    # answer = small_parsimony_unrooted(n, adjacency_list_str)
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

    answer = small_parsimony_unrooted(n, adjacency_list_str)

    with open(current_dir / "Wk3_2_output.txt", "w") as output_file:
        output_file.write(formatterer(answer))

