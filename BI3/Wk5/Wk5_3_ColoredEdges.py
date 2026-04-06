#########################################################################################
# Genome Sequencing  - Bioinformatic III Course from Coursera
#
# Week 4 - Genomes to the Breakpoint Graph
# Code Challenge: Implement ColoredEdges.
# Input: A genome P.
# Output: The collection of colored edges in the genome graph of P in the form (x, y).
# This is labeling COLORED edges, so just label the chromo2cycle as we did with negatives
# and just find the colored edges nodes :) 
"""
ColoredEdges(P)
     Edges ← an empty set
     for each chromosome Chromosome in P
          Nodes ← ChromosomeToCycle(Chromosome)
          for j ← 1 to |Chromosome|
               add the edge (Nodes2j, Nodes2j +1) to Edges
     return Edges
"""
# NOTE: YOU HAVE TO RUN THIS FROM THE BI3 FOLDER WITH THE COMMAND:
# python -m Wk5.Wk5_3_ColoredEdges (NO .py)
#########################################################################################

# Main code

import sys
from pathlib import Path as partho
sys.path.insert(0, str(partho(__file__).parent))
from Wk5.Wk5_1_Chromo2Cycle import chromosome_to_cycle


def colored_edges(P):
    """
    Finds all colored edges in the genome graph of P.
    For each chromosome in P, applies ChromosomeToCycle and then adds edges
    between adjacent cycle nodes (Nodes_{2j}, Nodes_{2j+1}) for j = 1 to |Chromosome|,
    wrapping around so the last node connects back to the first.

    Parameters:
    P -> genome string with one or more chromosomes e.g. "(+1 -2 -3)(+4 +5 -6)"

    Returns:
    list of tuples representing colored edges e.g. [(2, 4), (3, 6), (5, 1), ...]
    """
    chromosomes = ['(' + c for c in P.split('(') if c]
    # print(f"{chromosomes=}")
    edges = []
    for chromosome in chromosomes:
        nodes = chromosome_to_cycle(chromosome)
        # print("My nodes are: ", nodes)
        n = len(nodes)
        for j in range(1, n // 2 + 1):
            edges.append((nodes[2 * j - 1], nodes[(2 * j) % n]))
    return edges


###########################################################################

if __name__ == "__main__":
    # Sample test
    P_raw = "(+1 -2 -3)(+4 +5 -6)"
    # Expected answer = (2, 4), (3, 6), (5, 1), (8, 9), (10, 12), (11, 7)
    answer = colored_edges(P_raw)
    print(", ".join(str(e) for e in answer))

    # # From file
    # # Get dataset
    # from pathlib import Path as partho

    # current_dir = partho(__file__).parent
    # filename = input("Please enter the filename: ")
    # file_path = current_dir / filename

    # with open(file_path, "r") as file:
    #     P_raw = file.read().strip()

    # answer = colored_edges(P_raw)
    # # print(answer)
    # # answer = "(" + " ".join(map(str, answer)) + ")"
    # # answer = "(" + " ".join(f"{n:+d}" for n in answer) + ")"
    # answer = ", ".join(str(e) for e in answer)

    # with open("Wk4\Wk5_3_output.txt", "w") as output_file:
    #     output_file.write(str(answer))


    # Exam
    # permutation = ["+6", "-12", "-9", "+17", "+18", "-4", "+5",
    #                "-3", "+11", "+19", "+20", "+10", "+8", "+15", "-14",
    #                 "-13", "+2", "+7", "-16", "-1"]
    # # Expected answer = 8
    # answer = number_of_breakpoints(permutation)
    # print(answer)
