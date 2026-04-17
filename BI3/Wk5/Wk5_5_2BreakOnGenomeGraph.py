#########################################################################################
# Genome Sequencing  - Bioinformatic III Course from Coursera
#
# Week 4 - Genomes to the Breakpoint Graph
# Code Challenge: Implement GraphToGenome.
# Input: The colored edges ColoredEdges of a genome graph.
# Output: The genome P corresponding to this genome graph.
"""
2-BreakOnGenomeGraph(GenomeGraph, i1 , i2 , i3 , i4)
     remove colored edges (i1, i2) and (i3, i4) from GenomeGraph
     add colored edges (i1, i3) and (i2, i4) to GenomeGraph
     return GenomeGraph
"""
#########################################################################################

# Main code

import re


def parse_genome_graph(genome_graph_raw):
    """
    Parse a genome graph string into a list of edge tuples.

    Parameters:
        genome_graph_raw (str): Comma-separated tuples, e.g. "(2, 4), (3, 8)".

    Returns:
        list of tuple: List of (int, int) edge pairs.
    """
    pairs = re.findall(r'\((\d+),\s*(\d+)\)', genome_graph_raw)
    # print("Pairs are bhaiyya:", pairs)
    return [(int(a), int(b)) for a, b in pairs]


def two_break_on_genome_graph(genome_graph_raw, i1, i2, i3, i4):
    """
    Perform a 2-break operation on a genome graph.

    Removes colored edges (i1, i2) and (i3, i4) from the genome graph and
    adds colored edges (i1, i3) and (i2, i4) in their place.

    Parameters:
        genome_graph_raw (str): String representation of the genome graph's colored edges.
        i1 (int): First node of the first edge to remove.
        i2 (int): Second node of the first edge to remove.
        i3 (int): First node of the second edge to remove.
        i4 (int): Second node of the second edge to remove.

    Returns:
        str: String representation of the modified genome graph after the 2-break.
    """
    edges = parse_genome_graph(genome_graph_raw)
    # print("Returned edges yaar: ", edges)

    # Remove colored edges (i1, i2) and (i3, i4)
    edges = [e for e in edges if set(e) not in ({i1, i2}, {i3, i4})]
    # print("Edges with colored edges removed: ", edges)

    # Add colored edges (i1, i3) and (i2, i4)
    edges += [(i1, i3), (i2, i4)]
    # print("Edges with colored edges added: ", edges)

    return edges


###########################################################################

if __name__ == "__main__":
    # # Sample test
    # genome_graph_raw = "(2, 4), (3, 8), (7, 5), (6, 1)"
    # i1, i2, i3, i4 = 1, 6, 3, 8
    # # Expected answer = (2, 4), (3, 1), (7, 5), (6, 8)
    # answer = two_break_on_genome_graph(genome_graph_raw, i1, i2, i3, i4)
    # print(", ".join(f"({a}, {b})" for a, b in answer))

    # From file
    # Get dataset
    from pathlib import Path as partho

    current_dir = partho(__file__).parent
    filename = input("Please enter the filename: ")
    file_path = current_dir / filename

    with open(file_path, "r") as file:
        genome_graph_raw = file.readline().strip()
        i1, i2, i3, i4 = map(int, file.readline().strip().split(","))

    answer = two_break_on_genome_graph(genome_graph_raw, i1, i2, i3, i4)
    # print(answer)
    # answer = "(" + " ".join(map(str, answer)) + ")"
    # answer = "(" + " ".join(f"{n:+d}" for n in answer) + ")"
    # answer = ", ".join(str(e) for e in answer)
    # answer = "".join("(" + " ".join(f"{n:+d}" for n in chrom) + ")" for chrom in answer)
    print(", ".join(f"({a}, {b})" for a, b in answer))

    with open(current_dir / "Wk5_5_output.txt", "w") as output_file:
        answer = ", ".join(f"({a}, {b})" for a, b in answer)
        output_file.write(str(answer))

