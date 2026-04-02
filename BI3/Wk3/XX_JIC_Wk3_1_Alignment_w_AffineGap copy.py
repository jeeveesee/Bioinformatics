#########################################################################################
# Genome Sequencing  - Bioinformatic III Course from Coursera
#
# Week 3 - Alignment with Affine Gap Penalties Problem
# Code Challenge: Solve the Alignment with Affine Gap Penalties Problem.
# Input: A match reward, a mismatch penalty, a gap opening penalty, a gap extension penalty,
# and two nucleotide strings.
# Output: The maximum alignment score between v and w,
# followed by an alignment of v and w achieving this maximum score.

# NOTE: YOU HAVE TO RUN THIS FROM THE BI3 FOLDER WITH THE COMMAND:
# python -m Wk3.Wk3_1_Alignment_w_AffineGap (NO .py)
#########################################################################################
"""
lower(i, j) = max
                {
                    lower(i-1, j)  -epsilon <- gap extension
                    middle(i-1, j)  -sigma <- gap opening
                }
upper(i, j) = max
                {
                    upper(i, j-1)  -epsilon <- gap extension
                    middle(i, j-1)  -sigma <- gap opening
                }
middle(i, j) = max
                {
                    lower(i, j)
                    middle(i-1, j-1)  + score(vi, wj)
                    upper(i, j)
                }
"""
#########################################################################################

def affine_alignment(reward, mismatch, gap_opening, gap_ext, v, w):
    """
    Global alignment with affine gap penalties.

    Parameters:
        reward: Match reward value
        mismatch: Mismatch penalty
        gap_opening: Gap opening penalty
        gap_ext: Gap extension penalty
        v, w: Two sequences to align

    Returns:
        (score, aligned_v, aligned_w): Maximum alignment score and aligned sequences
    """
    n, m = len(v), len(w)

    # State names for clarity
    MATCH = 'MATCH'                      # Match/mismatch (diagonal move)
    GAP_IN_W_LOWER = 'GAP_IN_W_LOWER'    # Gap in w (vertical move, consumes v)
    GAP_IN_V_UPPER = 'GAP_IN_V_UPPER'    # Gap in v (horizontal move, consumes w)

    # Initialize matrices
    # gap_in_w_lower = vertical gaps (consuming from v, inserting gap in w)
    # gap_in_v_upper = horizontal gaps (consuming from w, inserting gap in v)
    # match = match/mismatch (diagonal)
    gap_in_w_lower = [[float('-inf')] * (m + 1) for _ in range(n + 1)]
    gap_in_v_upper = [[float('-inf')] * (m + 1) for _ in range(n + 1)]
    match = [[0] * (m + 1) for _ in range(n + 1)]

    # Traceback pointers
    # 0 = extended existing gap, 1 = opened new gap
    gap_in_w_lower_ptr = [[None] * (m + 1) for _ in range(n + 1)]
    gap_in_v_upper_ptr = [[None] * (m + 1) for _ in range(n + 1)]
    match_ptr = [[None] * (m + 1) for _ in range(n + 1)]

    # Initialize boundaries
    gap_in_w_lower[0][0] = gap_in_v_upper[0][0] = 0
    for i in range(1, n + 1):
        gap_in_w_lower[i][0] = -(gap_opening + (i - 1) * gap_ext)
        match[i][0] = gap_in_w_lower[i][0]
        gap_in_v_upper[i][0] = -(gap_opening + (i - 1) * gap_ext)

    for j in range(1, m + 1):
        gap_in_v_upper[0][j] = -(gap_opening + (j - 1) * gap_ext)
        match[0][j] = gap_in_v_upper[0][j]
        gap_in_w_lower[0][j] = -(gap_opening + (j - 1) * gap_ext)

    # Fill matrices
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # Gap in w (vertical move down - consume from v, add gap to w)
            extend_gap = gap_in_w_lower[i - 1][j] - gap_ext
            open_gap = match[i - 1][j] - gap_opening
            if extend_gap >= open_gap:
                gap_in_w_lower[i][j], gap_in_w_lower_ptr[i][j] = extend_gap, 0
            else:
                gap_in_w_lower[i][j], gap_in_w_lower_ptr[i][j] = open_gap, 1

            # Gap in v (horizontal move right - consume from w, add gap to v)
            extend_gap = gap_in_v_upper[i][j - 1] - gap_ext
            open_gap = match[i][j - 1] - gap_opening
            if extend_gap >= open_gap:
                gap_in_v_upper[i][j], gap_in_v_upper_ptr[i][j] = extend_gap, 0
            else:
                gap_in_v_upper[i][j], gap_in_v_upper_ptr[i][j] = open_gap, 1

            # Match/mismatch (diagonal move)
            diag_score = match[i - 1][j - 1] + (reward if v[i - 1] == w[j - 1] else -mismatch)

            # Find best option
            best_score = gap_in_w_lower[i][j]
            best_state = GAP_IN_W_LOWER
            if diag_score > best_score:
                best_score = diag_score
                best_state = MATCH
            if gap_in_v_upper[i][j] > best_score:
                best_score = gap_in_v_upper[i][j]
                best_state = GAP_IN_V_UPPER

            match[i][j] = best_score
            match_ptr[i][j] = best_state

    # Backtrack
    def handle_match_state(i, j):
        """Handle MATCH state - diagonal move consuming both sequences."""
        if i > 0 and j > 0:
            return v[i-1], w[j-1], i-1, j-1, match_ptr[i-1][j-1]
        elif i > 0:
            return v[i-1], '-', i-1, j, GAP_IN_W_LOWER
        else:
            return '-', w[j-1], i, j-1, GAP_IN_V_UPPER
    
    def handle_gap_in_w_lower_state(i, j):
        """Handle GAP_IN_W_LOWER state - vertical move consuming v, adding gap to w."""
        if i > 0:
            was_opened = gap_in_w_lower_ptr[i][j] == 1
            new_state = match_ptr[i-1][j] if was_opened else GAP_IN_W_LOWER
            return v[i-1], '-', i-1, j, new_state
        else:
            return '-', w[j-1], i, j-1, GAP_IN_V_UPPER
    
    def handle_gap_in_v_upper_state(i, j):
        """Handle GAP_IN_V_UPPER state - horizontal move consuming w, adding gap to v."""
        if j > 0:
            was_opened = gap_in_v_upper_ptr[i][j] == 1
            new_state = match_ptr[i][j-1] if was_opened else GAP_IN_V_UPPER
            return '-', w[j-1], i, j-1, new_state
        else:
            return v[i-1], '-', i-1, j, GAP_IN_W_LOWER
    
    def backtrack():
        alignment_1, alignment_2 = [], []
        i, j = n, m
        state = match_ptr[i][j] or MATCH
        
        while i > 0 or j > 0:
            if state == MATCH:
                char_v, char_w, i, j, state = handle_match_state(i, j)
            elif state == GAP_IN_W_LOWER:
                char_v, char_w, i, j, state = handle_gap_in_w_lower_state(i, j)
            elif state == GAP_IN_V_UPPER:
                char_v, char_w, i, j, state = handle_gap_in_v_upper_state(i, j)
            else:
                break
            
            alignment_1.append(char_v)
            alignment_2.append(char_w)
        
        return ''.join(reversed(alignment_1)), ''.join(reversed(alignment_2))
    
    aligned_v, aligned_w = backtrack()
    return match[n][m], aligned_v, aligned_w



###########################################################################

if __name__ == "__main__":

# Sample test
    reward = 1
    mismatch = 3
    gap_opening = 2
    gap_ext = 1
    v = "GA"
    w = "GTTA"
    # Expected answer =
    # -1
    # G--A
    # GTTA
    alignment_score, alignment_1, alignment_2 = affine_alignment(reward, mismatch, gap_opening, gap_ext, v, w)
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