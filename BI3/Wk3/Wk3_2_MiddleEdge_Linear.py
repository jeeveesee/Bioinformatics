#########################################################################################
# Genome Sequencing  - Bioinformatic III Course from Coursera
#
# Week 3 - Space efficient sequence alignment - Middle Edge in Linear Space
# Code Challenge: Solve the Middle Edge in Linear Space Problem (for protein strings).
# Input: A match reward, a mismatch penalty, an indel penalty, and two nucleotide strings.
# Output: A middle edge in the alignment graph.

# NOTE: YOU MAY HAVE TO RUN THIS FROM THE BI3 FOLDER WITH THE COMMAND:
# python -m Wk3.Wk3_2_MiddleEdge_Linear (NO .py)
#########################################################################################

def middle_edge(reward, mismatch, indel, v, w):
    # Get lengths of v and w
    n = len(v)
    m = len(w)

    # Get middle column index
    middle_col_index = m // 2

    # Initialize score arrays
    score_from_source = [0] * (n + 1)
    score_to_sink = [0] * (n + 1)

    # Fill score_from_source using forward dynamic programming
    for i in range(1, n + 1):
        score_from_source[i] = score_from_source[i - 1] - indel

    for j in range(1, middle_col_index + 1):
        prev_score_from_source = score_from_source[:]
        score_from_source[0] -= indel
        for i in range(1, n + 1):
            match_mismatch_score = reward if v[i - 1] == w[j - 1] else -mismatch
            score_from_source[i] = max(
                prev_score_from_source[i] - indel,      # Deletion
                score_from_source[i - 1] - indel,      # Insertion
                prev_score_from_source[i - 1] + match_mismatch_score  # Match/Mismatch
            )

    # Fill score_to_sink using backward dynamic programming
    for i in range(n - 1, -1, -1):
        score_to_sink[i] = score_to_sink[i + 1] - indel

    # We'll compute backward scores for each column and capture the arrays
    # for column middle_col_index (k) and column k+1 so we can decide
    # which outgoing edge (down/right/diagonal) from (i,k) lies on
    # an optimal path.
    score_to_sink_k = None       # scores for column k
    score_to_sink_kplus1 = None  # scores for column k+1

    for j in range(m - 1, -1, -1):
        prev_score_to_sink = score_to_sink[:]
        score_to_sink[n] -= indel
        for i in range(n - 1, -1, -1):
            match_mismatch_score = reward if v[i] == w[j] else -mismatch
            score_to_sink[i] = max(
                prev_score_to_sink[i] - indel,      # Deletion
                score_to_sink[i + 1] - indel,      # Insertion
                prev_score_to_sink[i + 1] + match_mismatch_score  # Match/Mismatch
            )
        # After updating for column j, score_to_sink corresponds to column j
        if j == middle_col_index + 1:
            score_to_sink_kplus1 = score_to_sink[:]
        if j == middle_col_index:
            score_to_sink_k = score_to_sink[:]

    # If middle_col_index is the last column, we may not have set k+1; use
    # the base scores (which represent column m) in that case.
    if score_to_sink_kplus1 is None:
        score_to_sink_kplus1 = [0] * (n + 1)
    if score_to_sink_k is None:
        score_to_sink_k = score_to_sink[:]

    # Find the maximum combined score and corresponding edge
    max_score = float('-inf')
    middle_edge_coord_1 = None
    middle_edge_coord_2 = None

    for i in range(n + 1):
        # Combined best score going through node (i, k)
        combined_score = score_from_source[i] + score_to_sink_k[i]
        if combined_score > max_score:
            max_score = combined_score
            middle_edge_coord_1 = (i, middle_col_index)
            # Decide which outgoing edge from (i,k) is on an optimal path.
            # Compute three candidate totals for right, down and diagonal.
            # Right: (i, k) -> (i, k+1)
            right_total = score_from_source[i] - indel + score_to_sink_kplus1[i]
            # Down: (i, k) -> (i+1, k)
            if i < n:
                down_total = score_from_source[i] - indel + score_to_sink_k[i + 1]
            else:
                down_total = float('-inf')
            # Diagonal: (i, k) -> (i+1, k+1)
            if i < n and middle_col_index < m:
                match_mismatch_score = reward if v[i] == w[middle_col_index] else -mismatch
                diag_total = score_from_source[i] + match_mismatch_score + score_to_sink_kplus1[i + 1]
            else:
                diag_total = float('-inf')

            # Choose the best move. Break ties in this priority:
            # right (horizontal) > down (vertical) > diagonal.
            if right_total >= down_total and right_total >= diag_total:
                middle_edge_coord_2 = (i, middle_col_index + 1)
            elif down_total >= diag_total:
                middle_edge_coord_2 = (i + 1, middle_col_index)
            else:
                middle_edge_coord_2 = (i + 1, middle_col_index + 1)

    return middle_edge_coord_1, middle_edge_coord_2

###########################################################################

if __name__ == "__main__":

# Sample test
    reward = 1
    mismatch = 1
    indel = 2
    v = "GAT"
    w = "GAGA"
    # Expected answer =
    # 2 2
    # 2 3
    middle_edge_coord_1, middle_edge_coord_2 = middle_edge(reward, mismatch, indel, v, w)
    # print('\n'.join(answer))
    print(*middle_edge_coord_1)
    print(*middle_edge_coord_2)

# # From file

#     # Get dataset
#     from pathlib import Path as partho
#     current_dir = partho(__file__).parent
#     filename = input("Please enter the filename: ")
#     file_path = current_dir / filename

#     #NOTE: For Wk3_2, there is a newline character after the peptide string, so remove that
#     with open(file_path, 'r') as file:
#         # Read in n and m
#         reward, mismatch, indel = file.readline().split()
#         w = file.readline().strip()
#         v = file.readline().strip()
#         # Down = [list(map(int, line.split())) for line in cleaned_input[1:separator]]

#     middle_edge_coord_1, middle_edge_coord_2 = middle_edge(int(reward), int(mismatch), int(indel), v, w)
#     # print('\n'.join(answer))
#     # print(answer)

#     with open("Wk3_2_output.txt", "w") as output_file:
#         output_file.write(" ".join(map(str, middle_edge_coord_1)))
#         output_file.write("\n")
#         output_file.write(" ".join(map(str, middle_edge_coord_2)))

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