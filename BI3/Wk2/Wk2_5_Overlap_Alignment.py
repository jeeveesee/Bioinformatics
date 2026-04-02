#########################################################################################
# Genome Sequencing  - Bioinformatic III Course from Coursera
#
# Week 2 - Overlap Alignment
# Code Challenge: Solve the Overlap Alignment Problem.
# Input: A match reward, a mismatch penalty, an indel penalty, and two nucleotide strings v and w.
# Output: The score of an optimal overlap alignment of v and w, followed by an alignment of a suffix
# v' of v and a prefix w' of w achieving this maximum score.

# NOTE: YOU HAVE TO RUN THIS FROM THE BI3 FOLDER WITH THE COMMAND:
# python -m Wk2.Wk2_2_Local_Alignment (NO .py)
#########################################################################################

def overlap_alignment(v, w, match_reward, mismatch_penalty, indel_sigma):
    """
    Uses overlap alignment to find the alignment between a suffix of string v and a prefix of string w by applying given match rewards, mismatch penalties, and indel penalties

    Parameters:
    match_reward -> reward for a match
    mismatch_penalty -> penalty for a mismatch
    indel_sigma -> indel penalty
    v, w -> two sequences

    Returns:
    Maximum score, alignment sequence 1, and alignment sequence 2 along with indels and mismatches
    """
    # Initialize scoring matrix
    n = len(v)
    m = len(w)
    S = [[0 for j in range(m + 1)] for i in range(n + 1)]

    # Initialize first row
    for j in range(m + 1):
        S[0][j] = 0

    # Initialize first column (all zeros for overlap alignment - can start at any suffix of v)
    for i in range(1, n + 1):
        S[i][0] = 0

    # Fill in the scoring matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if v[i - 1] == w[j - 1]:
                cost = match_reward
            else:
                cost = -mismatch_penalty
            S[i][j] = max(S[i - 1][j] - indel_sigma,      # Deletion
                          S[i][j - 1] - indel_sigma,      # Insertion
                          S[i - 1][j - 1] + cost)         # Match/Mismatch

    # Find max score in last row
    max_score = float('-inf')
    max_j = 0
    for j in range(m + 1):
        if S[n][j] > max_score:
            max_score = S[n][j]
            max_j = j

    # Traceback to get alignment
    alignment_1 = ""
    alignment_2 = ""
    i = n
    j = max_j

    while i > 0 and j > 0:
        if S[i][j] == S[i - 1][j] - indel_sigma:
            alignment_1 = v[i - 1] + alignment_1
            alignment_2 = "-" + alignment_2
            i -= 1
        elif S[i][j] == S[i][j - 1] - indel_sigma:
            alignment_1 = "-" + alignment_1
            alignment_2 = w[j - 1] + alignment_2
            j -= 1
        else:
            alignment_1 = v[i - 1] + alignment_1
            alignment_2 = w[j - 1] + alignment_2
            i -= 1
            j -= 1

    return max_score, alignment_1, alignment_2

###########################################################################

if __name__ == "__main__":

    # # Sample test
    # match_reward = 1
    # mismatch_penalty = 1
    # indel_sigma = 2

    # v = "GAGA"
    # w = "GAT"
    # # Expected answer =
    # # 2
    # # GA
    # # GA
    # score, alignment_1, alignment_2 = overlap_alignment(v, w, match_reward, mismatch_penalty, indel_sigma)
    # # print('\n'.join(answer))
    # print(score)
    # print(alignment_1)
    # print(alignment_2)

# # From file

#     # Get dataset
#     from pathlib import Path as partho
#     current_dir = partho(__file__).parent
#     filename = input("Please enter the filename: ")
#     file_path = current_dir / filename

#     #NOTE: For Wk3_2, there is a newline character after the peptide string, so remove that
#     with open(file_path, 'r') as file:
#         # Read in n and m
#         match_reward, mismatch_penalty, indel_sigma = file.readline().split()
#         v = file.readline().strip()
#         w = file.readline().strip()
#         # Down = [list(map(int, line.split())) for line in cleaned_input[1:separator]]

#     alignment_score, alignment_1, alignment_2 = overlap_alignment(v, w, int(match_reward), int(mismatch_penalty), int(indel_sigma))
#     # print('\n'.join(answer))
#     # print(answer)

#     with open("Wk2_5_output.txt", "w") as output_file:
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
    match_reward = 1
    mismatch_penalty = 1
    indel_sigma = 2

    v = "GAGA" #"AGTACATCAGAGGAGTT-ACATACTAACG"
    w = "GAT" #"AGTTCACAGGCTA-CGTACAGATATTACGACAGGCAGA"
    # Expected answer =
    # 2
    # GA
    # GA
    score, alignment_1, alignment_2 = overlap_alignment(v, w, match_reward, mismatch_penalty, indel_sigma)
    # print('\n'.join(answer))
    print(score)
    print(alignment_1)
    print(alignment_2)