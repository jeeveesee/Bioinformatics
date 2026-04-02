#########################################################################################
# Genome Sequencing  - Bioinformatic III Course from Coursera
#
# Week 2 - Global Alignment problem
# Code Challenge: Solve the Global Alignment Problem.
# Input: A match reward, a mismatch penalty, an indel penalty, and two nucleotide strings.
# Output: The maximum alignment score of these strings followed by
# an alignment achieving this maximum score.
# NOTE: YOU HAVE TO RUN THIS FROM THE BI3 FOLDER WITH THE COMMAND:
# python -m Wk2.Wk2_1_Global_Alignment (NO .py)
#########################################################################################

"""
s(i,j) = max(s(i-1, j) ---> -sigma
             s(i, j-1) ---> -sigma
             s(i-1, j-1) ---> +1, if v(i) = w(j)
             s(i-1, j-1) ---> -mu, if v(i) != w(j)
            )
"""
#########################################################################################

def global_alignment(v, w, reward, mismatch, indel_sigma):
    """
    Uses global alignment to find the alignment between any two strings by applying rewards, mismatch and indel penalties

    Parameters:
    reward -> reward value
    mismatch -> mismatch penalty
    indel_sigma -> indel penalty
    v, w -> two sequences

    Returns:
    Maximum score, alignment sequence 1, and alignment sequence 2 along with indels and mismatches
    """
    # Initialize scoring matrix
    n = len(v)
    m = len(w)
    S = [[0 for j in range(m + 1)] for i in range(n + 1)]

    # Initialize first row and column
    for i in range(1, n + 1):
        S[i][0] = S[i - 1][0] - indel_sigma
    for j in range(1, m + 1):
        S[0][j] = S[0][j - 1] - indel_sigma

    # Fill in the scoring matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = S[i - 1][j - 1] + (reward if v[i - 1] == w[j - 1] else -mismatch)
            delete = S[i - 1][j] - indel_sigma
            insert = S[i][j - 1] - indel_sigma
            S[i][j] = max(match, delete, insert)

    # Backtrack to find the optimal alignment
    alignment_1 = []
    alignment_2 = []
    i, j = n, m

    while i > 0 or j > 0:
        current_score = S[i][j]
        # Go diagonally
        if i > 0 and j > 0 and (current_score == S[i - 1][j - 1] + (reward if v[i - 1] == w[j - 1] else -mismatch)):
            alignment_1.append(v[i - 1])
            alignment_2.append(w[j - 1])
            i -= 1
            j -= 1
        # Go up (i changes, j remains the same)
        elif i > 0 and (current_score == S[i - 1][j] - indel_sigma):
            alignment_1.append(v[i - 1])
            alignment_2.append('-')
            i -= 1
        # Go left (i remains the same, j changes)
        else:
            alignment_1.append('-')
            alignment_2.append(w[j - 1])
            j -= 1

    alignment_1.reverse()
    alignment_2.reverse()

    return S[n][m], ''.join(alignment_1), ''.join(alignment_2)

###########################################################################

if __name__ == "__main__":

# # Sample test
#     reward = 1
#     mismatch = 1
#     indel_sigma = 2
#     v = "GAGA"
#     w = "GAT"
#     # Expected answer =
#     # -1
#     # GAGA
#     # GA-T
#     alignment_score, alignment_1, alignment_2 = global_alignment(v, w, reward, mismatch, indel_sigma)
#     # print('\n'.join(answer))
#     print(alignment_score)
#     print(alignment_1)
#     print(alignment_2)

# From file

    # Get dataset
    from pathlib import Path as partho
    current_dir = partho(__file__).parent
    filename = input("Please enter the filename: ")
    file_path = current_dir / filename

    #NOTE: For Wk3_2, there is a newline character after the peptide string, so remove that
    with open(file_path, 'r') as file:
        # Read in n and m
        reward, mismatch, indel_sigma = file.readline().split()
        v = file.readline().strip()
        w = file.readline().strip()
        # Down = [list(map(int, line.split())) for line in cleaned_input[1:separator]]

    alignment_score, alignment_1, alignment_2 = global_alignment(v, w, int(reward), int(mismatch), int(indel_sigma))
    # print('\n'.join(answer))
    # print(answer)

    with open("Wk2_1_output.txt", "w") as output_file:
        output_file.write(str(alignment_score))
        output_file.write("\n")
        output_file.write(str(alignment_1))
        output_file.write("\n")
        output_file.write(str(alignment_2))


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