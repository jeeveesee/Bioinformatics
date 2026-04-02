#########################################################################################
# Genome Sequencing  - Bioinformatic III Course from Coursera
#
# Week 1 - Longest Common Subsequence problem (OutputLCS)
# Code Challenge: Use OutputLCS (reproduced below) to solve the Longest Common Subsequence Problem.
# Input: Two strings s and t.
# Output: A longest common subsequence of s and t.
# (Note: more than one solution may exist, in which case you may output any one.)
#########################################################################################

"""
# 1 - Recursive

LCSBackTrack(v, w)
    for i ← 0 to |v|
        si, 0 ← 0
    for j ← 0 to |w|
        s0, j ← 0
    for i ← 1 to |v|
        for j ← 1 to |w|
            match ← 0
            if vi-1 = wj-1
                match ← 1
            si, j ← max{si-1, j , si,j-1 , si-1, j-1 + match }
            if si,j = si-1,j
                Backtracki, j ← "↓"
            else if si, j = si, j-1
                Backtracki, j ← "→"
            else if si, j = si-1, j-1 + match
                Backtracki, j ← "↘"
    return Backtrack

OutputLCS(backtrack, v, i, j)
    if i = 0 or j = 0
        return ""
    if backtracki, j = "↓"
        return OutputLCS(backtrack, v, i - 1, j)
    else if backtracki, j = "→"
        return OutputLCS(backtrack, v, i, j - 1)
    else
        return OutputLCS(backtrack, v, i - 1, j - 1) + vi
"""

"""
# 2 - Iterative
IterativeOutputLCS(Backtrack, v, w)
   LCS ← an empty string
    i ← length of string
    j ← length of string w
    while i > 0 and j > 0
        if Backtrack(i, j) = "↓"
            i ← i-1
        else if Backtrack(i,j) = "→"
            j ← j-1
        else if Backtrack(i,j) = "↘"
            LCS ← concatenate v[i] with LCS
            i ← i-1
            j ← j-1
    return LCS
"""
#########################################################################################

def output_LCS_iterative(v, w):
    """
    Uses dynamic programming to find the value of the longest common subsequence

    Parameters:
    s, t -> two peptide strings
    Down, Right -> Down and Right matrix values

    Returns:
    Length of longest path from source (0, 0) to sink (n, m)
    """
    LCS = ""
    i = len(v)
    j = len(w)
    Backtrack = LCS_Backtrack(v, w)

    while i > 0 and j > 0:
        if Backtrack[i][j] == "↓":
            i -= 1
        elif Backtrack[i][j] == "→":
            j -= 1
        elif Backtrack[i][j] == "↘":
            LCS = v[i - 1] + LCS
            i -= 1
            j -= 1
    return LCS


def LCS_Backtrack(v, w):
    """
    Finds the backtracking edges vector to in the longest common subpeptide problem

    Parameters:
    v, w -> two peptide strings

    Returns:
    Backtracked edges vector
    """
    n = len(v)
    m = len(w)
    S = [[0 for j in range(m + 1)] for i in range(n + 1)]
    # print(f"{S=}")
    Backtrack = [['' for j in range(m + 1)] for i in range(n + 1)]
    # print(f"{Backtrack=}")
    # for i in range(n + 1):
    #     S[i][0] = 0
    # for j in range(m + 1):
    #     S[0][j] = 0
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = 0
            if v[i - 1] == w[j - 1]:
                match = 1
            # si, j ← max{si-1, j , si,j-1 , si-1, j-1 + match }
            scores = [S[i - 1][j], S[i][j - 1], S[i - 1][j - 1] + match]
            S[i][j] = max(scores)
            if S[i][j] == S[i - 1][j]:
                Backtrack[i][j] = "↓"
            elif S[i][j] == S[i][j - 1]:
                Backtrack[i][j] = "→"
            elif S[i][j] == S[i - 1][j - 1] + match:
                Backtrack[i][j] = "↘"
    # print(Backtrack)
    return Backtrack

###########################################################################

if __name__ == "__main__":

# # Sample test
#     v = "AACCTTGG"
#     w = "ACACTGTGA"
#     # Expected answer = AACTGG
#     # backtrack = LCS_Backtrack(v, w)
#     answer = output_LCS_iterative(v, w)
#     # print('\n'.join(answer))
#     print(answer)

# Exam
    v = "AGACTG"
    w = "GTACGA"
    # Expected answer = AACTGG
    # backtrack = LCS_Backtrack(v, w)
    answer = output_LCS_iterative(v, w)
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
#         # Read in n and m
#         v = file.readline()
#         w = file.readline()
#         # Down = [list(map(int, line.split())) for line in cleaned_input[1:separator]]

#     answer = output_LCS_iterative(v, w)
#     # print('\n'.join(answer))
#     # print(answer)

#     with open("Wk1_3_output.txt", "w") as output_file:
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