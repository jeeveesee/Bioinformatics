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

    # STEP 1: The graph will store the adjacency list and the graph items, the nodes
    # Having a default dict is clutch! This is so that if a key is not present, we don't
    # get a keyerror; instead it makes an empty list for the key e.g.,
    # if 4 is not a key and you end up on 4, you get graph['4'] = []
    # instead of a key error
    graph = defaultdict(list)
    for line in adj_list:
        node, edges = line.split(": ")
        graph[node] = edges.split(" ")

    print(f"Nodes: {node}")
    print(f"Edges: {edges}")
    print(f"graph: {graph}")

    # STEP 2: Count in-degree and out-degree for each node
    indeg = Counter()
    outdeg = Counter()
    total_edges = 0

    for from_node, to_node in graph.items():
        outdeg[from_node] += len(to_node)
        total_edges += len(to_node)
        for v in to_node:
            indeg[v] += 1

    print(f"Out Degrees: {outdeg}")
    print(f"In Degrees: {indeg}")
    print(f"Total Edges: {total_edges}")

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

    # if no start node is found above (it's a circuit), pick any random outgoing node
    # if start_node is None:
    #     if outdeg:
    #         start_node = next(iter(outdeg.keys()))
    #     else:
    #         return []

    # Find keys with in-degrees > 1+ out-degree to make it the last node
    # keys_with_out_greater = []
    # keys_with_OUT_greater = [key for key in outdeg.keys() if outdeg[key] == indeg[key] + 1]
    # keys_with_IN_greater = [key for key in indeg.keys() if indeg[key] == outdeg[key] + 1]

    # for key in outdeg.keys():
    #     # Check if the out-degree value is exactly 1 greater than the in-degree value
    #     if outdeg[key] == indeg[key] + 1:
    #         keys_with_out_greater.append(key)

    # print(f"Keys with OUT greater: {keys_with_OUT_greater}")
    # print(f"Keys with IN greater: {keys_with_IN_greater}")

    # if len(keys_with_out_greater) > 1:
    #     raise ValueError("Too many options for starting key.")
    # elif len(keys_with_out_greater) == 0:
    #     print("No eligible starting key, choosing randomly")
    #     keys_with_out_greater = outdeg.key()[0]
    # else:
    #     keys_with_out_greater = keys_with_out_greater

    # Now let's run this through an Eulerian Cycle

    # Initialize variables for Hierholzer's algorithm
    # Cycle provides the Eulerian cycle
    # Stack provides a place where unexplored edges within the cycle exist
    # and start the next cycle at this point

    current_node = start_node or next(iter(outdeg.keys()))  # Start from an arbitrary node
    print(f"Starting Node: {current_node}")
    # print(f"Is it an integer?: {current_node == 0}")
    # print(f"Is it an alpha?: {current_node == '0'}")

    cycle = []
    stack = []

    while stack or graph[current_node]:
        print("\n")
        print(f"CURRENT CURRENT NODE: {current_node}")
        print(f"Graph[current_node]: {graph[current_node]}")
        print(f"NOT Graph[current_node]: {not graph[current_node]}")
        if not graph[current_node]:
            cycle.append(current_node)
            print(f"Cycle is: {cycle}")
            current_node = stack.pop()
            print(f"Current Node: {current_node}")
        else:
            stack.append(current_node)
            print(f"Stack: {stack}")
            next_node = graph[current_node].pop()
            print(f"Next Node, whcih will be current node: {next_node}")
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

#     answer = eulerian_cycle(adj_list.strip().split("\n"))
#     # print(answer)

#     with open("Wk2_1_output.txt", "w") as output_file:
#         output_file.write(answer)