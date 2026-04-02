#########################################################################################
# Genome Sequencing  - Bioinformatic II Course from Coursera
#
# Week 2 - Eulerian Cycles
# Input: The adjacency list of an Eulerian directed graph.
# Output: An Eulerian cycle in this graph.
#########################################################################################

# Function to read the adjacency list and contruct an Eulerian cycle
def eulerian_cycle(adj_list):
    # Create a dictionary to store the adjacency list
    graph = {}
    for line in adj_list:
        node, edges = line.split(": ")
        graph[node] = edges.split(" ")

    # print(f"Nodes: {node}")
    # print(f"Edges: {edges}")
    print(f"graph: {graph}")

    # Initialize variables for Hierholzer's algorithm
    # Cycle provides the Eulerian cycle
    # Stack provides a place where unexplored edges within the cycle exist
    # and start the next cycle at this point
    cycle = []
    stack = []
    current_node = list(graph.keys())[0]  # Start from an arbitrary node
    # print(f"Starting Node: {current_node}")
    # print(f"Is it an integer?: {current_node == 0}")
    # print(f"Is it an alpha?: {current_node == '0'}")

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
            next_node = graph[current_node].pop()
            # print(f"Next Node, whcih will be current node: {next_node}")
            current_node = next_node

    cycle.append(current_node)
    cycle.reverse()  # Reverse to get the correct order

    return "->".join(cycle)

if __name__ == "__main__":

# Sample test
#     adj_list = """
# 0: 3
# 1: 0
# 2: 1 6
# 3: 2
# 4: 2
# 5: 4
# 6: 5 8
# 7: 9
# 8: 7
# 9: 6
# """


    # answer = eulerian_cycle(adj_list.strip().split("\n"))
    # print(answer)

# From file

    # Get dataset
    from pathlib import Path as partho
    current_dir = partho(__file__).parent
    filename = input("Please enter the filename: ")
    file_path = current_dir / filename

    with open(file_path, 'r') as file:
       adj_list = file.read()

    answer = eulerian_cycle(adj_list.strip().split("\n"))
    # print(answer)

    with open("Wk2_1_output.txt", "w") as output_file:
        output_file.write(answer)