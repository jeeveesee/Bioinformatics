#########################################################################################
# Genome Sequencing  - Bioinformatic II Course from Coursera
#
# Week 2 - Eulerian Path
# Input: The adjacency list of a directed graph that has an Eulerian path
# Output: An Eulerian path in this graph.
# NOTE: This is the PATH, NOT a CYCLE!!!!
#########################################################################################

# Imports
from collections import Counter, defaultdict

# Function to read the adjacency list and contruct an Eulerian path
def eulerian_path(adj_list):
    """
    Creates an Eulerian path from an adjacency list
    :param adj_list -> adjacency list in a dictionary format
    :return -> Eulerian path
    """

    # STEP 1: The graph will store the adjacency list and the graph items, the nodes
    # Having a default dict is clutch! This is so that if a key is not present, we don't
    # get a keyerror; instead it makes an empty list for the key e.g.,
    # if 4 is not a key and you end up on 4, you get graph['4'] = []
    # instead of a key error
    graph = defaultdict(list)
    for line in adj_list:
        node, edges = line.split(": ")
        graph[node] = edges.split(" ")

    # print(f"Nodes: {node}")
    # print(f"Edges: {edges}")
    # print(f"graph: {graph}")

    # STEP 2: Count in-degree and out-degree for each node
    indeg = Counter()
    outdeg = Counter()
    total_edges = 0

    for from_node, to_node in graph.items():
        outdeg[from_node] += len(to_node)
        total_edges += len(to_node)
        for v in to_node:
            indeg[v] += 1

    # print(f"Out Degrees: {outdeg}")
    # print(f"In Degrees: {indeg}")
    # print(f"Total Edges: {total_edges}")

    if total_edges == 0:
        raise ValueError("No edges: a nontrivial Eulerian cycle does not exist.")

    # STEP 3: Find the starting node:
    # a. Either the node where out-degree = 1 + out_degree, OR
    # b. any random node if none available
    start_node = None
    for node in graph.keys():
        if outdeg[node] == indeg[node] + 1:
            start_node = node
            break

    current_node = start_node or next(iter(outdeg.keys()))  # Start from an arbitrary node
    # print(f"Starting Node: {current_node}")
    # print(f"Is it an integer?: {current_node == 0}")
    # print(f"Is it an alpha?: {current_node == '0'}")

    cycle = []
    stack = []

    while stack or graph[current_node]:
        # print("\n")
        # print(f"CURRENT CURRENT NODE: {current_node}")
        # print(f"Graph[current_node]: {graph[current_node]}")
        # print(f"NOT Graph[current_node]: {not graph[current_node]}")
        if not graph[current_node]:
            cycle.append(current_node)
            # print(f"Cycle is: {cycle}")
            current_node = stack.pop()
            # print(f"Current Node: {current_node}")
        else:
            stack.append(current_node)
            # print(f"Stack: {stack}")
            # Use FIFO order for edges so we follow the adjacency list order
            # (graph[current_node] is a list of neighbors in input order).
            # Pop from the front rather than the back to preserve that order.
            next_node = graph[current_node].pop(0)
            # print(f"Next Node, which will be current node: {next_node}")
            current_node = next_node

    cycle.append(current_node)
    cycle.reverse()  # Reverse to get the correct order

    return "->".join(cycle)

if __name__ == "__main__":

# Eulerian path with dead-ends
    adj_list = """
0: 2
1: 3
2: 1
3: 0 4
6: 3 7
7: 8
8: 9
9: 6
"""
    champaran = adj_list.strip().split("\n")
    print(champaran)
    answer = eulerian_path(champaran)
    print(answer)

# # From file

#     # Get dataset
#     from pathlib import Path as partho
#     current_dir = partho(__file__).parent
#     filename = input("Please enter the filename: ")
#     file_path = current_dir / filename

#     with open(file_path, 'r') as file:
#        adj_list = file.read()

#     answer = eulerian_path(adj_list.strip().split("\n"))
#     print(answer)

#     with open("Wk2_2_output.txt", "w") as output_file:
#         output_file.write(answer)