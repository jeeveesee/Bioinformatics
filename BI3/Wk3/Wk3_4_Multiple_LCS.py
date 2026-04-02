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
    lx, ly, lz = len(x), len(y), len(z)

    # Build 3D DP table
    # dp[i][j][k] = LCS score for x[:i], y[:j], z[:k]
    dp = [[[0] * (lz + 1) for _ in range(ly + 1)] for _ in range(lx + 1)]

    for i in range(1, lx + 1):
        for j in range(1, ly + 1):
            for k in range(1, lz + 1):
                # All 7 transitions
                candidates = [
                    dp[i-1][j][k],       # advance x only (gap in y, z)
                    dp[i][j-1][k],       # advance y only
                    dp[i][j][k-1],       # advance z only
                    dp[i-1][j-1][k],     # advance x, y
                    dp[i-1][j][k-1],     # advance x, z
                    dp[i][j-1][k-1],     # advance y, z
                    # advance all
                    dp[i-1][j-1][k-1] + (1 if x[i-1] == y[j-1] == z[k-1] else 0),
                ]
                dp[i][j][k] = max(candidates)

    score = dp[lx][ly][lz]

    # Backtrack to find alignment
    ax, ay, az = [], [], []
    i, j, k = lx, ly, lz

    while i > 0 or j > 0 or k > 0:
        # Determine which transition led to current cell
        cur = dp[i][j][k]

        if (i > 0 and j > 0 and k > 0
                and cur == (dp[i-1][j-1][k-1]
                            + (1 if x[i-1] == y[j-1] == z[k-1] else 0))):
            ax.append(x[i-1])
            ay.append(y[j-1])
            az.append(z[k-1])
            i -= 1
            j -= 1
            k -= 1
        elif i > 0 and j > 0 and cur == dp[i-1][j-1][k]:
            ax.append(x[i-1])
            ay.append(y[j-1])
            az.append('-')
            i -= 1
            j -= 1
        elif i > 0 and k > 0 and cur == dp[i-1][j][k-1]:
            ax.append(x[i-1])
            ay.append('-')
            az.append(z[k-1])
            i -= 1
            k -= 1
        elif j > 0 and k > 0 and cur == dp[i][j-1][k-1]:
            ax.append('-')
            ay.append(y[j-1])
            az.append(z[k-1])
            j -= 1
            k -= 1
        elif i > 0 and cur == dp[i-1][j][k]:
            ax.append(x[i-1])
            ay.append('-')
            az.append('-')
            i -= 1
        elif j > 0 and cur == dp[i][j-1][k]:
            ax.append('-')
            ay.append(y[j-1])
            az.append('-')
            j -= 1
        else:
            ax.append('-')
            ay.append('-')
            az.append(z[k-1])
            k -= 1

    # Reverse since we backtracked
    ax = ''.join(reversed(ax))
    ay = ''.join(reversed(ay))
    az = ''.join(reversed(az))

    return score, ax, ay, az


###########################################################################

if __name__ == "__main__":
    # # Sample test
    # x = "ATATCCG"
    # y = "TCCGA"
    # z = "ATGTACTG"
    # # Expected answer =
    # # Score: 3
    # # ATATCC-G-
    # # ---TCC-GA
    # # ATGTACTG-
    # alignment_score, alignment_1, alignment_2, alignment_3 = multiple_lcs(x, y, z)
    # print(alignment_score)
    # print(alignment_1)
    # print(alignment_2)
    # print(alignment_3)

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

    # with open("Wk3_4_output.txt", "w") as output_file:
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

    # EXAM
    x = "TCTAGCGAAC"
    y = "ATTACCGATC"
    z = "TTCACTGACG"
    alignment_score, alignment_1, alignment_2, alignment_3 = multiple_lcs(x, y, z)
    print(alignment_score)
    print(alignment_1)
    print(alignment_2)
    print(alignment_3)