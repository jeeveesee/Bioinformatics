#########################################################################################
# Genome Sequencing  - Bioinformatic III Course from Coursera
#
# Week 3 - Alignment with Affine Gap Penalties Problem
# Code Challenge: Solve the Alignment with Affine Gap Penalties Problem.
# Input: A match reward, a mismatch penalty, a gap opening penalty, a gap extension penalty,
# and two nucleotide strings.
# Output: The maximum alignment score between v and w,
# followed by an alignment of v and w achieving this maximum score.
#########################################################################################
"""
Formulas from the lecture:

lower(i, j) = max {
    lower(i-1, j) - epsilon    <- extend gap in w
    middle(i-1, j) - sigma     <- open new gap in w
}

upper(i, j) = max {
    upper(i, j-1) - epsilon    <- extend gap in v
    middle(i, j-1) - sigma     <- open new gap in v
}

middle(i, j) = max {
    lower(i, j)                <- come from lower (gap in w)
    middle(i-1, j-1) + score(vi, wj)  <- match/mismatch
    upper(i, j)                <- come from upper (gap in v)
}
"""
#########################################################################################

def affine_alignment(reward, mismatch, gap_opening, gap_ext, v, w):
    """
    Global alignment with affine gap penalties.

    Parameters:
        reward: Match reward value
        mismatch: Mismatch penalty (positive value, will be subtracted)
        gap_opening: Gap opening penalty sigma (positive value, will be subtracted)
        gap_ext: Gap extension penalty epsilon (positive value, will be subtracted)
        v, w: Two sequences to align

    Returns:
        (score, aligned_v, aligned_w): Maximum alignment score and aligned sequences
    """
    n, m = len(v), len(w)
    epsilon = gap_ext
    sigma = gap_opening

    # Initialize matrices
    lower = [[float('-inf')] * (m + 1) for _ in range(n + 1)]
    upper = [[float('-inf')] * (m + 1) for _ in range(n + 1)]
    middle = [[0] * (m + 1) for _ in range(n + 1)]

    # Traceback pointers: which matrix did we come from?
    backtrack = [[None] * (m + 1) for _ in range(n + 1)]

    # Initialize first row and column
    lower[0][0] = upper[0][0] = 0
    for i in range(1, n + 1):
        lower[i][0] = -(sigma + (i - 1) * epsilon)
        middle[i][0] = lower[i][0]
        backtrack[i][0] = 'lower'

    for j in range(1, m + 1):
        upper[0][j] = -(sigma + (j - 1) * epsilon)
        middle[0][j] = upper[0][j]
        backtrack[0][j] = 'upper'

    # Fill matrices according to the formulas
    for i in range(1, n + 1):
        for j in range(1, m + 1):

            # LOWER(i,j) = max of two options:
            #   1. Extend existing gap: lower(i-1, j) - epsilon
            #   2. Open new gap: middle(i-1, j) - sigma
            # print("LOWER")
            # print("==========")
            # print("Extend gap value, ", lower[i-1][j] - epsilon)
            # print("Open gap value, ", middle[i-1][j] - sigma)
            lower[i][j] = max(
                lower[i-1][j] - epsilon,    # extend gap
                middle[i-1][j] - sigma       # open gap
            )

            # UPPER(i,j) = max of two options:
            #   1. Extend existing gap: upper(i, j-1) - epsilon
            #   2. Open new gap: middle(i, j-1) - sigma
            # print("UPPER")
            # print("==========")
            # print("Extend gap value, ", upper[i][j-1] - epsilon)
            # print("Open gap value, ", middle[i][j-1] - sigma)
            upper[i][j] = max(
                upper[i][j-1] - epsilon,     # extend gap
                middle[i][j-1] - sigma        # open gap
            )

            # MIDDLE(i,j) = max of three options:
            #   1. Come from lower (gap in w)
            #   2. Match/mismatch from diagonal
            #   3. Come from upper (gap in v)
            match_score = reward if v[i-1] == w[j-1] else -mismatch
            # print("DIAGONAL")
            # print("==========")
            # print(v[i-1] == w[j-1])
            # print(match_score)
            diagonal = middle[i-1][j-1] + match_score

            # Find which is best
            if lower[i][j] >= diagonal and lower[i][j] >= upper[i][j]:
                middle[i][j] = lower[i][j]
                backtrack[i][j] = 'lower'
            elif diagonal >= lower[i][j] and diagonal >= upper[i][j]:
                middle[i][j] = diagonal
                backtrack[i][j] = 'middle'
            else:
                middle[i][j] = upper[i][j]
                backtrack[i][j] = 'upper'

    # Backtrack to build alignment
    def backtrack_alignment():
        alignment_v, alignment_w = [], []
        i, j = n, m
        current = backtrack[i][j]

        while i > 0 or j > 0:
            if current == 'middle':
                # Came from diagonal - match or mismatch
                if i > 0 and j > 0:
                    alignment_v.append(v[i-1])
                    alignment_w.append(w[j-1])
                    i, j = i-1, j-1
                    current = backtrack[i][j] if i > 0 or j > 0 else None
                else:
                    # Edge case: consume remaining
                    if i > 0:
                        alignment_v.append(v[i-1])
                        alignment_w.append('-')
                        i -= 1
                    else:
                        alignment_v.append('-')
                        alignment_w.append(w[j-1])
                        j -= 1

            elif current == 'lower':
                # Gap in w (consume from v)
                if i > 0:
                    alignment_v.append(v[i-1])
                    alignment_w.append('-')
                    # Check if we stay in lower or go back to middle
                    if lower[i][j] == lower[i-1][j] - epsilon:
                        i -= 1
                        current = 'lower'  # extended gap
                    else:
                        i -= 1
                        current = backtrack[i][j] if i > 0 or j > 0 else None  # opened gap
                else:
                    alignment_v.append('-')
                    alignment_w.append(w[j-1])
                    j -= 1
                    current = 'upper'

            elif current == 'upper':
                # Gap in v (consume from w)
                if j > 0:
                    alignment_v.append('-')
                    alignment_w.append(w[j-1])
                    # Check if we stay in upper or go back to middle
                    if upper[i][j] == upper[i][j-1] - epsilon:
                        j -= 1
                        current = 'upper'  # extended gap
                    else:
                        j -= 1
                        current = backtrack[i][j] if i > 0 or j > 0 else None  # opened gap
                else:
                    alignment_v.append(v[i-1])
                    alignment_w.append('-')
                    i -= 1
                    current = 'lower'
            else:
                break

        return ''.join(reversed(alignment_v)), ''.join(reversed(alignment_w))

    aligned_v, aligned_w = backtrack_alignment()
    return middle[n][m], aligned_v, aligned_w


###########################################################################

if __name__ == "__main__":

# # Sample test
#     reward = 1
#     mismatch = 3
#     gap_opening = 2
#     gap_ext = 1
#     v = "GA"
#     w = "GTTA"
#     # Expected answer =
#     # -1
#     # G--A
#     # GTTA
#     alignment_score, alignment_1, alignment_2 = affine_alignment(reward, mismatch, gap_opening, gap_ext, v, w)
#     # print('\n'.join(answer))
#     print(alignment_score)
#     print(alignment_1)
#     print(alignment_2)

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

# EXAM
    reward = 1
    mismatch = 1
    gap_opening = 4
    gap_ext = 1
    v = "ACGTTAC"
    w = "ATGCAGT"
    alignment_score, alignment_1, alignment_2 = affine_alignment(reward, mismatch, gap_opening, gap_ext, v, w)
    # print('\n'.join(answer))
    print(alignment_score)
    print(alignment_1)
    print(alignment_2)