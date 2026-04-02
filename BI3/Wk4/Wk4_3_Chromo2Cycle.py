#########################################################################################
# Genome Sequencing  - Bioinformatic III Course from Coursera
#
# Week 4 - Genomes to the Breakpoint Graph
# Code Challenge: Implement ChromosomeToCycle.
# Input: A chromosome Chromosome containing n synteny blocks.
# Output: The sequence Nodes of integers between 1 and 2n
# resulting from applying ChromosomeToCycle to Chromosome.
# To do the naming, label all nodes as 1, 2, 3, 4, 5, 6 etc. then for negatives just switch
# the numbering... e.g., if the chromo is +1, -2, -3... label them 1, 2, 3, 4, 5, 6
# and then switch -2 ---> 4, 3 and -3 to 6, 5 for the final of 1, 2, 4, 3, 6, 5
"""
ChromosomeToCycle(Chromosome)
     for j ← 1 to |Chromosome|
          i ← Chromosomej
          if i > 0
               Nodes2j-1 ←2i-1
               Nodes2j ← 2i
          else
               Nodes2j-1 ← -2i
               Nodes2j ←-2i-1
     return Nodes
"""
#########################################################################################

# Main code


def chromosome_to_cycle(chromosome):
    """
    Converts a chromosome (list of signed integers) to its cycle node representation.
    For each synteny block i in the chromosome:
      - if i > 0: assigns nodes 2i-1 and 2i
      - if i < 0: assigns nodes -2i and -2i-1

    Parameters:
    chromosome -> list of signed integers e.g. [1, -2, -3, 4]

    Returns:
    list of integers representing the cycle nodes e.g. [1, 2, 4, 3, 6, 5, 7, 8]
    """
    chromosome = [int(x) for x in chromosome.strip("()").split()]
    nodes = [0] * (2 * len(chromosome))
    for j, i in enumerate(chromosome, start=1):
        if i > 0:
            nodes[2 * j - 2] = 2 * i - 1  # Nodes_{2j-1}
            nodes[2 * j - 1] = 2 * i      # Nodes_{2j}
        else:
            nodes[2 * j - 2] = -2 * i     # Nodes_{2j-1}
            nodes[2 * j - 1] = -2 * i - 1 # Nodes_{2j}
    return nodes


###########################################################################

if __name__ == "__main__":
    # # Sample test
    # chromosome_raw = "(+1 -2 -3 +4)"
    # # Expected answer = (1 2 4 3 6 5 7 8)
    # answer = chromosome_to_cycle(chromosome_raw)
    # print("(" + " ".join(map(str, answer)) + ")")

    # From file

    # Get dataset
    from pathlib import Path as partho

    current_dir = partho(__file__).parent
    filename = input("Please enter the filename: ")
    file_path = current_dir / filename

    with open(file_path, "r") as file:
        chromosome = file.read().strip()

    answer = chromosome_to_cycle(chromosome)
    # print(answer)
    answer = "(" + " ".join(map(str, answer)) + ")"

    with open("Wk4_3_output.txt", "w") as output_file:
        output_file.write(str(answer))

    # Exam
    # permutation = ["+6", "-12", "-9", "+17", "+18", "-4", "+5",
    #                "-3", "+11", "+19", "+20", "+10", "+8", "+15", "-14",
    #                 "-13", "+2", "+7", "-16", "-1"]
    # # Expected answer = 8
    # answer = number_of_breakpoints(permutation)
    # print(answer)
