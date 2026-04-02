#########################################################################################
# Genome Sequencing  - Bioinformatic III Course from Coursera
#
# Week 2 - Local Alignment problem
# Code Challenge: Solve the Local Alignment Problem.
# Input: Two protein strings written in the single-letter amino acid alphabet.
# Output: The maximum score of a local alignment of the strings,
# followed by a local alignment of these strings achieving the maximum score.
# Use the PAM250 scoring matrix for matches and mismatches
# as well as the indel penalty σ = 5.
# NOTE: YOU HAVE TO RUN THIS FROM THE BI3 FOLDER WITH THE COMMAND:
# python -m Wk2.Wk2_2_Local_Alignment (NO .py)
#########################################################################################

"""
s(i,j) = max(0
             s(i-1, j) ---> -sigma
             s(i, j-1) ---> -sigma
             s(i-1, j-1) ---> Pam250 score
            )
"""
#########################################################################################

PAM250 = {'A': {'A': 2, 'C': -2, 'E': 0, 'D': 0, 'G': 1, 'F': -3, 'I': -1, 'H': -1, 'K': -1, 'M': -1, 'L': -2, 'N': 0, 'Q': 0, 'P': 1, 'S': 1, 'R': -2, 'T': 1, 'W': -6, 'V': 0, 'Y': -3},
          'C': {'A': -2, 'C': 12, 'E': -5, 'D': -5, 'G': -3, 'F': -4, 'I': -2, 'H': -3, 'K': -5, 'M': -5, 'L': -6, 'N': -4, 'Q': -5, 'P': -3, 'S': 0, 'R': -4, 'T': -2, 'W': -8, 'V': -2, 'Y': 0},
          'E': {'A': 0, 'C': -5, 'E': 4, 'D': 3, 'G': 0, 'F': -5, 'I': -2, 'H': 1, 'K': 0, 'M': -2, 'L': -3, 'N': 1, 'Q': 2, 'P': -1, 'S': 0, 'R': -1, 'T': 0, 'W': -7, 'V': -2, 'Y': -4},
          'D': {'A': 0, 'C': -5, 'E': 3, 'D': 4, 'G': 1, 'F': -6, 'I': -2, 'H': 1, 'K': 0, 'M': -3, 'L': -4, 'N': 2, 'Q': 2, 'P': -1, 'S': 0, 'R': -1, 'T': 0, 'W': -7, 'V': -2, 'Y': -4},
          'G': {'A': 1, 'C': -3, 'E': 0, 'D': 1, 'G': 5, 'F': -5, 'I': -3, 'H': -2, 'K': -2, 'M': -3, 'L': -4, 'N': 0, 'Q': -1, 'P': 0, 'S': 1, 'R': -3, 'T': 0, 'W': -7, 'V': -1, 'Y': -5},
          'F': {'A': -3, 'C': -4, 'E': -5, 'D': -6, 'G': -5, 'F': 9, 'I': 1, 'H': -2, 'K': -5, 'M': 0, 'L': 2, 'N': -3, 'Q': -5, 'P': -5, 'S': -3, 'R': -4, 'T': -3, 'W': 0, 'V': -1, 'Y': 7},
          'I': {'A': -1, 'C': -2, 'E': -2, 'D': -2, 'G': -3, 'F': 1, 'I': 5, 'H': -2, 'K': -2, 'M': 2, 'L': 2, 'N': -2, 'Q': -2, 'P': -2, 'S': -1, 'R': -2, 'T': 0, 'W': -5, 'V': 4, 'Y': -1},
          'H': {'A': -1, 'C': -3, 'E': 1, 'D': 1, 'G': -2, 'F': -2, 'I': -2, 'H': 6, 'K': 0, 'M': -2, 'L': -2, 'N': 2, 'Q': 3, 'P': 0, 'S': -1, 'R': 2, 'T': -1, 'W': -3, 'V': -2, 'Y': 0},
          'K': {'A': -1, 'C': -5, 'E': 0, 'D': 0, 'G': -2, 'F': -5, 'I': -2, 'H': 0, 'K': 5, 'M': 0, 'L': -3, 'N': 1, 'Q': 1, 'P': -1, 'S': 0, 'R': 3, 'T': 0, 'W': -3, 'V': -2, 'Y': -4},
          'M': {'A': -1, 'C': -5, 'E': -2, 'D': -3, 'G': -3, 'F': 0, 'I': 2, 'H': -2, 'K': 0, 'M': 6, 'L': 4, 'N': -2, 'Q': -1, 'P': -2, 'S': -2, 'R': 0, 'T': -1, 'W': -4, 'V': 2, 'Y': -2},
          'L': {'A': -2, 'C': -6, 'E': -3, 'D': -4, 'G': -4, 'F': 2, 'I': 2, 'H': -2, 'K': -3, 'M': 4, 'L': 6, 'N': -3, 'Q': -2, 'P': -3, 'S': -3, 'R': -3, 'T': -2, 'W': -2, 'V': 2, 'Y': -1},
          'N': {'A': 0, 'C': -4, 'E': 1, 'D': 2, 'G': 0, 'F': -3, 'I': -2, 'H': 2, 'K': 1, 'M': -2, 'L': -3, 'N': 2, 'Q': 1, 'P': 0, 'S': 1, 'R': 0, 'T': 0, 'W': -4, 'V': -2, 'Y': -2},
          'Q': {'A': 0, 'C': -5, 'E': 2, 'D': 2, 'G': -1, 'F': -5, 'I': -2, 'H': 3, 'K': 1, 'M': -1, 'L': -2, 'N': 1, 'Q': 4, 'P': 0, 'S': -1, 'R': 1, 'T': -1, 'W': -5, 'V': -2, 'Y': -4},
          'P': {'A': 1, 'C': -3, 'E': -1, 'D': -1, 'G': 0, 'F': -5, 'I': -2, 'H': 0, 'K': -1, 'M': -2, 'L': -3, 'N': 0, 'Q': 0, 'P': 6, 'S': 1, 'R': 0, 'T': 0, 'W': -6, 'V': -1, 'Y': -5},
          'S': {'A': 1, 'C': 0, 'E': 0, 'D': 0, 'G': 1, 'F': -3, 'I': -1, 'H': -1, 'K': 0, 'M': -2, 'L': -3, 'N': 1, 'Q': -1, 'P': 1, 'S': 2, 'R': 0, 'T': 1, 'W': -2, 'V': -1, 'Y': -3},
          'R': {'A': -2, 'C': -4, 'E': -1, 'D': -1, 'G': -3, 'F': -4, 'I': -2, 'H': 2, 'K': 3, 'M': 0, 'L': -3, 'N': 0, 'Q': 1, 'P': 0, 'S': 0, 'R': 6, 'T': -1, 'W': 2, 'V': -2, 'Y': -4},
          'T': {'A': 1, 'C': -2, 'E': 0, 'D': 0, 'G': 0, 'F': -3, 'I': 0, 'H': -1, 'K': 0, 'M': -1, 'L': -2, 'N': 0, 'Q': -1, 'P': 0, 'S': 1, 'R': -1, 'T': 3, 'W': -5, 'V': 0, 'Y': -3},
          'W': {'A': -6, 'C': -8, 'E': -7, 'D': -7, 'G': -7, 'F': 0, 'I': -5, 'H': -3, 'K': -3, 'M': -4, 'L': -2, 'N': -4, 'Q': -5, 'P': -6, 'S': -2, 'R': 2, 'T': -5, 'W': 17, 'V': -6, 'Y': 0},
          'V': {'A': 0, 'C': -2, 'E': -2, 'D': -2, 'G': -1, 'F': -1, 'I': 4, 'H': -2, 'K': -2, 'M': 2, 'L': 2, 'N': -2, 'Q': -2, 'P': -1, 'S': -1, 'R': -2, 'T': 0, 'W': -6, 'V': 4, 'Y': -2},
          'Y': {'A': -3, 'C': 0, 'E': -4, 'D': -4, 'G': -5, 'F': 7, 'I': -1, 'H': 0, 'K': -4, 'M': -2, 'L': -1, 'N': -2, 'Q': -4, 'P': -5, 'S': -3, 'R': -4, 'T': -3, 'W': 0, 'V': -2, 'Y': 10}}

def local_alignment(v, w, indel_sigma):
    """
    Uses local alignment to find the alignment between any two strings by applying PAM250 scoring matrix and indel penalties

    Parameters:
    indel_sigma -> indel penalty
    v, w -> two sequences

    Returns:
    Maximum score, alignment sequence 1, and alignment sequence 2 along with indels and mismatches
    """
    # Initialize scoring matrix
    n = len(v)
    m = len(w)
    S = [[0 for j in range(m + 1)] for i in range(n + 1)]

    max_score = 0
    max_pos = (0, 0)

    # Fill in the scoring matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = S[i - 1][j - 1] + PAM250[v[i - 1]][w[j - 1]]
            delete = S[i - 1][j] - indel_sigma
            insert = S[i][j - 1] - indel_sigma
            S[i][j] = max(0, match, delete, insert)

            if S[i][j] > max_score:
                max_score = S[i][j]
                max_pos = (i, j)

    # Backtrack to find the optimal alignment
    alignment_1 = []
    alignment_2 = []
    i, j = max_pos

    while S[i][j] != 0:
        current_score = S[i][j]
        if current_score == S[i - 1][j] - indel_sigma: # Deletion
            alignment_1.append(v[i - 1])
            alignment_2.append('-')
            i -= 1
        elif current_score == S[i][j - 1] - indel_sigma: # Insertion
            alignment_1.append('-')
            alignment_2.append(w[j - 1])
            j -= 1
        else:
            alignment_1.append(v[i - 1]) # Mismatches and matches
            alignment_2.append(w[j - 1])
            i -= 1
            j -= 1

    alignment_1.reverse()
    alignment_2.reverse()

    return max_score, ''.join(alignment_1), ''.join(alignment_2)

###########################################################################

if __name__ == "__main__":

# Sample test
    v = "MEANLY"
    w = "PENALTY"
    indel_sigma = 5
    # Expected answer =
    # 15
    # EANL-Y
    # ENALTY
    alignment_score, alignment_1, alignment_2 = local_alignment(v, w, indel_sigma)
    # print('\n'.join(answer))
    print(alignment_score)
    print(alignment_1)
    print(alignment_2)

# # From file

#     # Get dataset
#     from pathlib import Path as partho
#     current_dir = partho(__file__).parent
#     filename = input("Please enter the filename: ")
#     file_path = current_dir / filename

#     #NOTE: For Wk3_2, there is a newline character after the peptide string, so remove that
#     with open(file_path, 'r') as file:
#         # Read in n and m
#         reward, mismatch, indel_sigma = file.readline().split()
#         v = file.readline().strip()
#         w = file.readline().strip()
#         # Down = [list(map(int, line.split())) for line in cleaned_input[1:separator]]

#     alignment_score, alignment_1, alignment_2 = global_alignment(v, w, int(reward), int(mismatch), int(indel_sigma))
#     # print('\n'.join(answer))
#     # print(answer)

#     with open("Wk2_1_output.txt", "w") as output_file:
#         output_file.write(str(alignment_score))
#         output_file.write("\n")
#         output_file.write(str(alignment_1))
#         output_file.write("\n")
#         output_file.write(str(alignment_2))


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