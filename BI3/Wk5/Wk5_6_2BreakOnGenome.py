#########################################################################################
# Genome Sequencing  - Bioinformatic III Course from Coursera
#
# Week 4 - Genomes to the Breakpoint Graph
# Code Challenge: Implement 2-BreakOnGenome.
# Input: A genome P, followed by indices i1 , i2 , i3 , and i4 .
# Output: The genome P' resulting from applying the 2-break operation
# 2-BreakOnGenome(GenomeGraph i1 , i2 , i3 , i4 ).
"""
2-BreakOnGenome(genome_raw, i1 , i2 , i3 , i4 )
     genome ← parse genome_raw
     GenomeGraph ← BlackEdges(genome) and ColoredEdges(genome)
     GenomeGraph ← 2-BreakOnGenomeGraph(GenomeGraph, i1 , i2 , i3 , i4 )
     P ← GraphToGenome(GenomeGraph)
     return P
"""
# NOTE: YOU HAVE TO RUN THIS FROM THE BI3 FOLDER WITH THE COMMAND:
# python -m Wk5.Wk5_6_2BreakOnGenome (NO .py)
#########################################################################################

# Main code

import sys
from pathlib import Path as partho

sys.path.insert(0, str(partho(__file__).parent))
from BI3.Wk5.Wk5_3_ColoredEdges import colored_edges
from BI3.Wk5.Wk5_4_Graph2Genome import graph_to_genome
from BI3.Wk5.Wk5_5_2BreakOnGenomeGraph import two_break_on_genome_graph


def two_break_on_genome(genome_raw, i1, i2, i3, i4):
    """
    Applies a 2-break operation to a genome.

    Converts the genome to its genome graph (colored edges), applies the
    2-break on the genome graph, then reconstructs the resulting genome P'.

    Parameters:
        genome_raw (str): Genome string, e.g. "(+1 -2 -4 +3)".
        i1 (int): First node of the first colored edge to remove.
        i2 (int): Second node of the first colored edge to remove.
        i3 (int): First node of the second colored edge to remove.
        i4 (int): Second node of the second colored edge to remove.

    Returns:
        list of tuple: Each tuple holds the signed-string elements of one
                       chromosome, e.g. [('+1', '-2'), ('-3', '+4')].
    """
    # Build genome graph from colored edges
    edges = colored_edges(genome_raw)

    # Apply 2-break on the genome graph
    edges_str = ", ".join(f"({a}, {b})" for a, b in edges)
    new_edges = two_break_on_genome_graph(edges_str, i1, i2, i3, i4)

    # Reconstruct genome from updated colored edges
    new_edges_str = ", ".join(f"({a}, {b})" for a, b in new_edges)
    chromosomes = graph_to_genome(new_edges_str)

    return chromosomes


###########################################################################

if __name__ == "__main__":
    # # Sample test
    # genome_raw = "(+1 -2 -4 +3)"
    # i1, i2, i3, i4 = 1, 6, 3, 8
    # # Expected answer = (+1 -2) (-3 +4)
    # answer = two_break_on_genome(genome_raw, i1, i2, i3, i4)
    # answer = [tuple(f"{n:+d}" for n in chrom) for chrom in answer]
    # print("".join(f"({a} {b})" for a, b in answer))

    # From file
    # Get dataset
    from pathlib import Path as partho

    current_dir = partho(__file__).parent
    filename = input("Please enter the filename: ")
    file_path = current_dir / filename

    with open(file_path, "r") as file:
        genome_raw = file.readline().strip()
        i1, i2, i3, i4 = map(int, file.readline().strip().split(","))

    answer = two_break_on_genome(genome_raw, i1, i2, i3, i4)
    answer = [tuple(f"{n:+d}" for n in chrom) for chrom in answer]
    answer = "".join("(" + " ".join(chrom) + ")" for chrom in answer)
    # print(answer)

    # answer = "(" + " ".join(map(str, answer)) + ")"
    # answer = "(" + " ".join(f"{n:+d}" for n in answer) + ")"
    # answer = ", ".join(str(e) for e in answer)
    # answer = "".join("(" + " ".join(f"{n:+d}" for n in chrom) + ")" for chrom in answer)

    with open(current_dir / "Wk4_8_output.txt", "w") as output_file:
        output_file.write(str(answer))


    # Exam
    # permutation = ["+6", "-12", "-9", "+17", "+18", "-4", "+5",
    #                "-3", "+11", "+19", "+20", "+10", "+8", "+15", "-14",
    #                 "-13", "+2", "+7", "-16", "-1"]
    # # Expected answer = 8
    # answer = number_of_breakpoints(permutation)
    # print(answer)
