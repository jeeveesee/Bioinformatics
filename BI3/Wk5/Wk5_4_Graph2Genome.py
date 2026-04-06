#########################################################################################
# Genome Sequencing  - Bioinformatic III Course from Coursera
#
# Week 4 - Genomes to the Breakpoint Graph
# Code Challenge: Implement GraphToGenome.
# Input: The colored edges ColoredEdges of a genome graph.
# Output: The genome P corresponding to this genome graph.
"""
GraphToGenome(GenomeGraph)
     P ← an empty set of chromosomes
     for each cycle Nodes in GenomeGraph
          Nodes ← sequence of nodes in this cycle (starting from node 1)
          Chromosome ← CycleToChromosome(Nodes)
          add Chromosome to P
     return P
"""
# NOTE: YOU HAVE TO RUN THIS FROM THE BI3 FOLDER WITH THE COMMAND:
# python -m Wk5.Wk5_4_Graph2Genome (NO .py)
#########################################################################################

# Main code

import sys
from pathlib import Path as partho
sys.path.insert(0, str(partho(__file__).parent))
from Wk5.Wk5_2_Cycle2Chromo import cycle_to_chromosome


def graph_to_genome(genome_graph_raw):
    """
    Converts the colored edges of a genome graph back to genome P.

    Parameters:
    genome_graph_raw -> str, colored edges e.g. "(2, 4), (3, 6), (5, 1), (7, 9), (10, 12), (11, 8)"

    Returns:
    list of chromosomes, each a list of signed integers e.g. [[1, -2, -3], [-4, 5, -6]]
    """
    colored_adj = {}
    for pair in genome_graph_raw.split("), ("):
        u, v = [int(x.strip("() ")) for x in pair.split(",")]
        colored_adj[u] = v
        colored_adj[v] = u

    def black(n):
        return n + 1 if n % 2 == 1 else n - 1

    p = []
    visited = set()

    for x in sorted(colored_adj):
        if x in visited:
            continue
        cycle_nodes = []
        current = x
        while True:
            cycle_nodes.append(current)
            visited.add(current)
            partner = black(current)
            cycle_nodes.append(partner)
            visited.add(partner)
            current = colored_adj[partner]
            if current == x:
                break
        p.append(cycle_to_chromosome("(" + " ".join(map(str, cycle_nodes)) + ")"))

    return p


###########################################################################

if __name__ == "__main__":
    # Sample test
    genome_graph_raw = "(2, 4), (3, 6), (5, 1), (7, 9), (10, 12), (11, 8)"
    # Expected answer = (+1 -2 -3)(-4 +5 -6)
    answer = graph_to_genome(genome_graph_raw)
    print("".join("(" + " ".join(f"{n:+d}" for n in chrom) + ")" for chrom in answer))

    # # From file
    # # Get dataset
    # from pathlib import Path as partho

    # current_dir = partho(__file__).parent
    # filename = input("Please enter the filename: ")
    # file_path = current_dir / filename

    # with open(file_path, "r") as file:
    #     genome_graph_raw = file.read().strip()

    # answer = graph_to_genome(genome_graph_raw)
    # # print(answer)
    # # answer = "(" + " ".join(map(str, answer)) + ")"
    # # answer = "(" + " ".join(f"{n:+d}" for n in answer) + ")"
    # # answer = ", ".join(str(e) for e in answer)
    # answer = "".join("(" + " ".join(f"{n:+d}" for n in chrom) + ")" for chrom in answer)

    # with open("Wk5\Wk5_4_output.txt", "w") as output_file:
    #     output_file.write(str(answer))


    # Exam
    # permutation = ["+6", "-12", "-9", "+17", "+18", "-4", "+5",
    #                "-3", "+11", "+19", "+20", "+10", "+8", "+15", "-14",
    #                 "-13", "+2", "+7", "-16", "-1"]
    # # Expected answer = 8
    # answer = number_of_breakpoints(permutation)
    # print(answer)
