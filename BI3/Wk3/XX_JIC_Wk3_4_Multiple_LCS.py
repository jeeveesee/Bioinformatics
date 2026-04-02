#########################################################################################
# Genome Sequencing  - Bioinformatic III Course from Coursera
#
# Week 3 - Multiple Longest Common Subsequence problem (LCS)
# The score of a column of the alignment matrix is equal to 1
# if all of the column's symbols are identical, and 0 if even one symbol disagrees
# Code Challenge: Solve the Multiple Longest Common Subsequence Problem.
# Input: Three DNA strings of length at most 10.
# Output: The length of a longest common subsequence of these three strings,
# followed by a multiple alignment of the three strings corresponding to such an alignment.

# NOTE: YOU MAY HAVE TO RUN THIS FROM THE BI3 FOLDER WITH THE COMMAND:
# python -m Wk3.Wk3_4_Multiple_LCS (NO .py)
#########################################################################################
"""
Max score:
==================
s(i, j, k) = max(
                    s(i-1, j, k) + score(x(i), -, -)
                    s(i, j-1, k) + score(-, y(j), -)
                    s(i, j, k-1) + score(-, -, z(k))
                    s(i-1, j-1, k) + score(x(i), y(j), -)
                    s(i-1, j, k-1) + score(x(i), -, z(k))
                    s(i, j-1, k-1) + score(-, y(j), z(k))
                    s(i-1, j-1, k-1) + score(x(i), y(j), z(k))
                )
"""
#########################################################################################

def multiple_lcs(x, y, z):
    """
    Returns the score of the longest common subsequence of three strings,
    and one such longest common subsequence.
    Args:
    Three DNA strings x, y, and z.
    Returns:
    The score of the longest common subsequence of x, y, and z,
    followed by a multiple alignment of x, y, and z corresponding to such an alignment.
    (If multiple solutions exist, you may return any one.)
    Sample Input:
    ATATCCG
    TCCGA
    ATGTACTG
    """
    # Create a 3D matrix to store the scores of the longest common subsequence of the three strings
    # Initialize the first layer of the matrix to 0
    # Fill in the rest of the matrix using dynamic programming
    # The score of a column of the alignment matrix is equal to 1 if all of the column's symbols are identical, and 0 if even one symbol disagrees
    # Backtrack through the matrix to find one of the longest common subsequences
    m = len(x)
    n = len(y)
    o = len(z)
    # Create a 3D matrix to store the scores of the longest common subsequence of the three strings
    score_matrix = [[[0 for _ in range(o + 1)] for _ in range(n + 1)] for _ in range(m + 1)]
    # Fill in the rest of the matrix using dynamic programming
    # Standard LCS recurrence for three strings (3-transition):
    # if all three chars equal -> dp[i][j][k] = dp[i-1][j-1][k-1] + 1
    # else -> max of dp[i-1][j][k], dp[i][j-1][k], dp[i][j][k-1]
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            for k in range(1, o + 1):
                if x[i - 1] == y[j - 1] == z[k - 1]:
                    score_matrix[i][j][k] = score_matrix[i - 1][j - 1][k - 1] + 1
                else:
                    score_matrix[i][j][k] = max(
                        score_matrix[i - 1][j][k],
                        score_matrix[i][j - 1][k],
                        score_matrix[i][j][k - 1],
                    )
    # Backtrack through the matrix to find one of the longest common subsequences
    alignment_1 = ""
    alignment_2 = ""
    alignment_3 = ""
    i = m
    j = n
    k = o
    # Simplified one-character-at-a-time backtracking (plus the full-match case)
    while i > 0 or j > 0 or k > 0:
        if i > 0 and j > 0 and k > 0 and x[i - 1] == y[j - 1] == z[k - 1]:
            alignment_1 = x[i - 1] + alignment_1
            alignment_2 = y[j - 1] + alignment_2
            alignment_3 = z[k - 1] + alignment_3
            i -= 1
            j -= 1
            k -= 1
        else:
            # Prefer consuming from x, then y, then z when scores tie
            if i > 0 and score_matrix[i][j][k] == score_matrix[i - 1][j][k]:
                alignment_1 = x[i - 1] + alignment_1
                alignment_2 = "-" + alignment_2
                alignment_3 = "-" + alignment_3
                i -= 1
            elif j > 0 and score_matrix[i][j][k] == score_matrix[i][j - 1][k]:
                alignment_1 = "-" + alignment_1
                alignment_2 = y[j - 1] + alignment_2
                alignment_3 = "-" + alignment_3
                j -= 1
            elif k > 0 and score_matrix[i][j][k] == score_matrix[i][j][k - 1]:
                alignment_1 = "-" + alignment_1
                alignment_2 = "-" + alignment_2
                alignment_3 = z[k - 1] + alignment_3
                k -= 1
            else:
                # Edge fallback: if only one index remains, consume it
                if i > 0 and j == 0 and k == 0:
                    alignment_1 = x[i - 1] + alignment_1
                    alignment_2 = "-" + alignment_2
                    alignment_3 = "-" + alignment_3
                    i -= 1
                elif j > 0 and i == 0 and k == 0:
                    alignment_1 = "-" + alignment_1
                    alignment_2 = y[j - 1] + alignment_2
                    alignment_3 = "-" + alignment_3
                    j -= 1
                elif k > 0 and i == 0 and j == 0:
                    alignment_1 = "-" + alignment_1
                    alignment_2 = "-" + alignment_2
                    alignment_3 = z[k - 1] + alignment_3
                    k -= 1
    # The score of the longest common subsequence is the value in the bottom right corner of the matrix
    alignment_score = score_matrix[m][n][o]
    return alignment_score, alignment_1, alignment_2, alignment_3


###########################################################################

if __name__ == "__main__":
    # Sample test
    x = "ATATCCG"
    y = "TCCGA"
    z = "ATGTACTG"
    # Expected answer =
    # Score: 3
    # ATATCC-G-
    # ---TCC-GA
    # ATGTACTG-
    alignment_score, alignment_1, alignment_2, alignment_3 = multiple_lcs(x, y, z)
    print(alignment_score)
    print(alignment_1)
    print(alignment_2)
    print(alignment_3)

    # # From file

    # # Get dataset
    # from pathlib import Path as partho
    # current_dir = partho(__file__).parent
    # filename = input("Please enter the filename: ")
    # file_path = current_dir / filename

    # #NOTE: For Wk3_2, there is a newline character after the peptide string, so remove that
    # with open(file_path, 'r') as file:
    #     # Read in x, y and z
    #     # reward, mismatch, indel = file.readline().split()
    #     x = file.readline().strip()
    #     y = file.readline().strip()
    #     z = file.readline().strip()
    #     # Down = [list(map(int, line.split())) for line in cleaned_input[1:separator]]

    # alignment_score, alignment_1, alignment_2, alignment_3 = multiple_lcs(x, y, z)
    # # print('\n'.join(answer))
    # # print(answer)

    # with open("Wk3/Wk3_4_output.txt", "w") as output_file:
    #     output_file.write(str(alignment_score))
    #     output_file.write("\n")
    #     output_file.write(alignment_1)
    #     output_file.write("\n")
    #     output_file.write(alignment_2)
    #     output_file.write("\n")
    #     output_file.write(alignment_3)
    #     # output_file.write(" ".join(map(str, middle_edge_coord_1)))
    #     # output_file.write("\n")

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