#########################################################################################
# Genome Sequencing  - Bioinformatic II Course from Coursera
#
# Week 2 - Maximal Non Branching Paths
#   Input: The adjacency list of a graph whose nodes are integers
#   Output: The collection of all maximal nonbranching paths in this graph
# Pseudocode:
"""
MaximalNonBranchingPaths(Graph)
    Paths ← empty list
    for each node v in Graph
        if v is not a 1-in-1-out node
            if out(v) > 0
                for each outgoing edge (v, w) from v
                    NonBranchingPath ← the path consisting of single edge (v, w)
                    while w is a 1-in-1-out node
                        extend NonBranchingPath by the edge (w, u) 
                        w ← u
                    add NonBranchingPath to the set Paths
    for each isolated cycle Cycle in Graph
        add Cycle to Paths
    return Paths
"""
#########################################################################################

def maximal_non_branching_paths(adj_list):
    """Return all maximal non-branching paths in the directed graph.

    :param adj: adjacency dict representing the directed graph
    :return: list of paths, each path is a list of nodes
    """
    from collections import defaultdict

    # Compute indegree and outdegree
    indeg = defaultdict(int)
    outdeg = defaultdict(int)
    adj_dict = defaultdict(list)

    # Turn adj_list into a dictionary if not already a dict
    if isinstance(adj_list, dict):
        adj_dict = adj_list
    else:
        for line in adj_list:
            node, edges = line.split(": ")
            adj_dict[node] = edges.split(" ")

    for src in adj_dict:
        outdeg[src] += len(adj_dict[src])
        for dst in adj_dict[src]:
            indeg[dst] += 1

    # print(f"{adj_dict=}")
    # print(f"{outdeg=}")
    # print(f"{indeg=}")
    paths = []

    # Function to check if a node is 1-in-1-out
    def is_1_in_1_out(node):
        return indeg[node] == 1 and outdeg[node] == 1

    # Find non-branching paths
    for v in adj_dict:
        # print(f"\n{v=}")
        if not is_1_in_1_out(v):
            if outdeg[v] > 0:
                for w in adj_dict[v]:
                    # print(f"{w=}")
                    non_branching_path = [v, w]
                    # print(f"{non_branching_path=}")
                    while is_1_in_1_out(w):
                        # print(f"{adj_dict[w]=}")
                        u = adj_dict[w][0] # Coz only the first one will be a non-branching one
                        # print(f"{u=}")
                        non_branching_path.append(u)
                        # print(f"{non_branching_path=}")
                        w = u
                    paths.append(non_branching_path)
                    # print(f"{paths=}")

    # Find isolated cycles
    # Only follow chains made entirely of 1-in-1-out nodes. This prevents
    # traversing into nodes that don't have outgoing edges and avoids
    # IndexError from indexing into empty adjacency lists.
    # print("----- Starting Isolated Cycles -----")
    visited = set()
    for v in adj_dict:
        # If it is not a 1-1 node, it means it is part of an isolated path
        if not is_1_in_1_out(v) or v in visited:
            continue
        # print(f"Isolated {v=}")
        cycle = [v]
        visited.add(v)
        # print(f"Isolated {cycle=}")
        # print(f"Isolated {visited=}")
        # start from the unique outgoing neighbor
        w = adj_dict[v][0]
        # print(f"Isolated {adj_dict=}")
        # print(f"Isolated {w=}")

        # Walk while the next node is also 1-in-1-out
        while True:
            # If next node is not 1-in-1-out, this isn't an isolated cycle
            if not is_1_in_1_out(w):
                break

            # If we came back to the start, close the cycle and add it
            if w == v:
                cycle.append(v)
                paths.append(cycle)
                break

            # If we've already visited this node (but it's not the start), abort
            if w in visited:
                break

            cycle.append(w)
            visited.add(w)
            # safe to index adj_dict[w][0] because is_1_in_1_out(w) is True
            w = adj_dict[w][0]

    return paths


if __name__ == "__main__":

# Sample dataset
    adj_list = """
1: 2
2: 3
3: 4 5
6: 7
7: 6
"""
    adj_list = adj_list.strip().split("\n")
    print(adj_list)
    answer = maximal_non_branching_paths(adj_list)
    print("\n", answer)

# # # From file
#     # Get dataset
#     from pathlib import Path as partho
#     current_dir = partho(__file__).parent
#     filename = input("Please enter the filename: ")
#     file_path = current_dir / filename

#     with open(file_path, 'r') as file:
#        adj_list = file.read()

#     chullu = maximal_non_branching_paths(int(k), int(d), read_pairs)
#     print(chullu)
#     # The above is done to close the loop and make it circular
#     # as the overlap is k-1, so you leave out k-1 nucleotides at the end


#     with open("Wk2_8_output.txt", "w") as output_file:
#         output_file.write(chullu)
