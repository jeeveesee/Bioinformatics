#########################################################################################
# Genome Sequencing  - Bioinformatic III Course from Coursera
#
# Week 1 - Longest Common Subsequence problem for ANY DAG
# Code Challenge: Solve the Longest Path in a DAG Problem
# Input: An integer representing the starting node to consider in a graph,
# followed by an integer representing the ending node to consider, followed
# by a list of edges in the graph. The edge notation "0 1 7" indicates that an edge connects node 0 to node 1 with weight 7.
# You may assume a given topological order corresponding to nodes in increasing order.
# Output: The length of a longest path in the graph, followed by a longest path as a sequence of space-separated node labels.
# (If multiple longest paths exist, you may return any one.)

# Generalized backtracking method -->
# Sb = max of all predecessors a of node b {Sa + weight of edge from a to b}

# (Note: more than one solution may exist, in which case you may output any one.)
#########################################################################################

def LCS_general(start_node: int, end_node: int, dag: list[list[int]]) -> tuple[float, list]:
    """
    Uses dynamic programming to find the value of the longest common subsequence
    between two strings v and w.

    Args:
        start_node: Integer representing the starting node in the DAG.
        end_node: Integer representing the ending node in the DAG.
        dag: List of edges in the DAG, where each edge is represented as a list [u, v, weight].
    Returns:
        A tuple containing the length of the longest path and the longest path as a list of node labels.
    """
    # Initialize distance and backtrack dictionaries
    distance = {node: float('-inf') for node in range(end_node + 1)}
    distance[start_node] = 0
    backtrack = {node: None for node in range(end_node + 1)}

    # Process each edge in the DAG
    for u, v, weight in dag:
        if distance[u] + weight > distance[v]:
            distance[v] = distance[u] + weight
            backtrack[v] = u

    # Reconstruct the longest path from start_node to end_node
    path = []
    current_node = end_node
    while current_node is not None:
        path.append(current_node)
        current_node = backtrack[current_node]
    path.reverse()

    return distance[end_node], path


###########################################################################

if __name__ == "__main__":

# # Sample test
#     start_node = 0
#     end_node = 4
#     dag = [[0, 1, 7], [0, 2, 4], [2, 3, 2], [1, 4, 1], [3, 4, 3]]
#     # Expected answer =
#     # 9
#     # 0 2 3 4
#     answer = LCS_general(start_node, end_node, dag)
#     # print('\n'.join(answer))
#     print(answer)

# Exam
    start_node = 1
    end_node = 7
    dag = [[1, 2, 5], [1, 3, 6], [1, 4, 5], [2, 3, 2],
           [2, 6, 9], [3, 5, 4], [3, 6, 3], [3, 7, 7],
           [4, 5, 4], [4, 6, 5], [5, 7, 2], [6, 7, 1]]
    # Expected answer =
    # 9
    # 0 2 3 4
    answer = LCS_general(start_node, end_node, dag)
    # print('\n'.join(answer))
    print(answer)

# # From file

#     # Get dataset
#     from pathlib import Path as partho
#     current_dir = partho(__file__).parent
#     filename = input("Please enter the filename: ")
#     file_path = current_dir / filename

#     #NOTE: For Wk3_2, there is a newline character after the peptide string, so remove that
#     with open(file_path, 'r') as file:
#         lines = file.readlines()
#         start_node, end_node = map(int, lines[0].split())
#         dag = [[int(x) for x in line.split()] for line in lines[1:]]
#         # Down = [list(map(int, line.split())) for line in cleaned_input[1:separator]]

#     answer = LCS_general(start_node, end_node, dag)
#     # print('\n'.join(answer))
#     # print(answer)
#     print(answer[0])
#     print(*answer[1])

#     with open("Wk1_4_output.txt", "w") as output_file:
#         output_file.write(str(answer))


# # For Exercise Break

#     # Get dataset
#     from pathlib import Path as partho
#     current_dir = partho(__file__).parent
#     filename = input("Please enter the filename: ")
#     file_path = current_dir / filename

#     #NOTE: For Wk3_2 exercise break, we need to concatenate all the kmers on different lines
#     with open(file_path, 'r') as file:
#        dna_string = file.read().replace('\n', '')
#     # print(dna_string)

#     # Val-Lys-Leu-Phe-Pro-Trp-Phe-Asn-Gln-Tyr
#     peptide_string = 'VKLFPWFNQY'

#     answer = peptide_encoding(RNA_CODON_MAP, dna_string, peptide_string)
#     # print(len(answer))

#     with open("Wk3_2_exercisebreak_output.txt", "w") as output_file:
#         output_file.write('\n'.join(answer))
#         # Answer is 0!!

# # EXAM
#     peptide_string = 'PEEP'
#     experimental_spectrum = [0, 97, 129, 129, 129, 194, 226, 323, 323, 355, 452]
#     answer = linear_scoring(peptide_string, experimental_spectrum)
#     # print('\n'.join(answer))
#     print(answer)