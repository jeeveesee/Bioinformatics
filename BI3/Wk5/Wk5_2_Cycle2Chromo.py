#########################################################################################
# Genome Sequencing  - Bioinformatic III Course from Coursera
#
# Week 4 - Genomes to the Breakpoint Graph
# Code Challenge: Implement CycleToChromosome.
# Input: A sequence Nodes of integers between 1 and 2n.
# Output: The chromosome Chromosome containing n synteny blocks
# resulting from applying CycleToChromosome to Nodes.
"""
CycleToChromosome(Nodes)
     for j ← 1 to |Nodes|/2
          if Nodes2j-1 < Nodes2j
               Chromosomej ← Nodes2j /2
          else
               Chromosomej ← -Nodes2j-1/2
     return Chromosome
"""
#########################################################################################

# Main code
def cycle_to_chromosome(nodes):
    """
    Converts a cycle node representation to its chromosome (list of signed integers).
    For each pair of nodes (Nodes_{2j-1}, Nodes_{2j}):
      - if Nodes_{2j-1} < Nodes_{2j}: assigns Chromosome_j = Nodes_{2j} / 2
      - else: assigns Chromosome_j = -Nodes_{2j-1} / 2

    Parameters:
    nodes -> list of signed integers e.g. (1 2 4 3 6 5 7 8)

    Returns:
    list of integers representing the circular chromosome e.g. [1, -2, -3, 4]

    """
    nodes = [int(x) for x in nodes.strip("()").split()]
    chromosome = []
    for j in range(1, len(nodes) // 2 + 1):
        if nodes[2 * j - 2] < nodes[2 * j - 1]:
            chromosome.append(nodes[2 * j - 1] // 2)
        else:
            chromosome.append(-nodes[2 * j - 2] // 2)
    return chromosome


###########################################################################

if __name__ == "__main__":
    # # Sample test
    # nodes_raw = "(1 2 4 3 6 5 7 8)"
    # # Expected answer = (+1 -2 -3 +4)
    # answer = cycle_to_chromosome(nodes_raw)
    # # print("(" + " ".join(map(str, answer)) + ")")
    # print("(" + " ".join(f"{n:+d}" for n in answer) + ")")
    # # print(answer)

    # From file
    # Get dataset
    from pathlib import Path as partho

    current_dir = partho(__file__).parent
    filename = input("Please enter the filename: ")
    file_path = current_dir / filename

    with open(file_path, "r") as file:
        nodes_raw = file.read().strip()

    answer = cycle_to_chromosome(nodes_raw)
    # print(answer)
    # answer = "(" + " ".join(map(str, answer)) + ")"
    answer = "(" + " ".join(f"{n:+d}" for n in answer) + ")"

    with open("Wk5_2_output.txt", "w") as output_file:
        output_file.write(str(answer))


    # Exam
    # permutation = ["+6", "-12", "-9", "+17", "+18", "-4", "+5",
    #                "-3", "+11", "+19", "+20", "+10", "+8", "+15", "-14",
    #                 "-13", "+2", "+7", "-16", "-1"]
    # # Expected answer = 8
    # answer = number_of_breakpoints(permutation)
    # print(answer)
