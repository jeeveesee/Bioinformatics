#########################################################################################
# Genome Sequencing  - Bioinformatic III Course from Coursera
#
# Week 1 - Dynamic Programming Manhattan Tourist
# Code Challenge: Find the length of a longest path in the Manhattan Tourist Problem.
# Input: Integers n and m, followed by an n × (m + 1) matrix Down and an (n + 1) × m matrix Right. The two matrices are separated by the "-" symbol.
# Output: The length of a longest path from source (0, 0) to sink (n, m) in the rectangular grid whose edges are defined by the matrices Down and Right.


"""
ManhattanTourist(n, m, Down, Right)
    s0, 0 ← 0
    for i ← 1 to n
        si, 0 ← si-1, 0 + downi-1, 0
    for j ← 1 to m
        s0, j ← s0, j−1 + right0, j-1
    for i ← 1 to n
        for j ← 1 to m
            si, j ← max{si - 1, j + downi-1, j, si, j - 1 + righti, j-1}
    return sn, m
"""
# NOTE: YOU HAVE TO RUN THIS FROM THE BI3 FOLDER WITH THE COMMAND:
# python -m Wk1.Wk1_2_ManhattanTourist (NO .py)
#########################################################################################

def Manhattan_Tourist(n, m, Down, Right):
    """
    Uses dynamic programming to find the value of the longest path for a manhattan tourist problem

    Parameters:
    n, m -> integers to make n x (m+1) Down matrix and (n+1) x m Right matrix
    Down, Right -> Down and Right matrix values

    Returns:
    Length of longest path from source (0, 0) to sink (n, m)
    """
    # s0, 0 ← 0
    S = [[0 for j in range(m + 1)] for i in range(n + 1)]
    # for i ← 1 to n
    for i in range(1, n + 1):
        # si, 0 ← si-1, 0 + downi-1, 0
        S[i][0] = S[i - 1][0] + Down[i - 1][0]
    # for j ← 1 to m
    for j in range(1, m + 1):
        # s0, j ← s0, j−1 + right0, j-1
        S[0][j] = S[0][j - 1] + Right[0][j - 1]
    # for i ← 1 to n
    for i in range(1, n + 1):
        # for j ← 1 to m
        for j in range(1, m + 1):
            # si, j ← max{si - 1, j + downi-1, j, si, j - 1 + righti, j-1}
            S[i][j] = max(S[i - 1][j] + Down[i - 1][j], S[i][j - 1] + Right[i][j - 1])
    # return sn, m
    return S[n][m]

###########################################################################

if __name__ == "__main__":

# Sample test
    # n = 4
    # m = 4
    # Down = [[1, 0, 2, 4, 3], [4, 6, 5, 2, 1], [4, 4, 5, 2, 1], [5, 6, 8, 5, 3]]
    # Right = [[3, 2, 4, 0], [3, 2, 4, 2], [0, 7, 3, 3], [3, 3, 0, 2], [1, 3, 2, 2]]
    # # Expected answer = 34
    # answer = Manhattan_Tourist(n, m, Down, Right)
    # # print('\n'.join(answer))
    # print(answer)

# From file

    # Get dataset
    from pathlib import Path as partho
    current_dir = partho(__file__).parent
    filename = input("Please enter the filename: ")
    file_path = current_dir / filename

    #NOTE: For Wk3_2, there is a newline character after the peptide string, so remove that
    with open(file_path, 'r') as file:
        # Read in n and m
        raw_input = file.readlines()
        cleaned_input = [item.strip() for item in raw_input]
        n, m = map(int, raw_input[0].split())
        # # Read in the two matrices
        separator = cleaned_input.index("-")
        Down = [list(map(int, line.split())) for line in cleaned_input[1:separator]]
        Right = [list(map(int, line.split())) for line in cleaned_input[separator + 1:]]

    answer = Manhattan_Tourist(n, m, Down, Right)
    # print('\n'.join(answer))
    # print(answer)

    with open("Wk1_2_output.txt", "w") as output_file:
        output_file.write(str(answer))


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