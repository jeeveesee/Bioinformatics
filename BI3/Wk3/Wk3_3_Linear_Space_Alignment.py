#########################################################################################
# Genome Sequencing  - Bioinformatic III Course from Coursera
#
# Week 3 - Linear space alignment
# Code Challenge: Implement LinearSpaceAlignment to solve the
# Global Alignment Problem for a large dataset.
# Input: A match reward, a mismatch penalty, an indel penalty,
# and two long (2000 nucleotide) DNA strings.
# Output: The maximum alignment score of these strings,
# followed by an alignment achieving this maximum score.


# NOTE: I wasn't able to run this and get the correct answer within Cogniterra
# The issue lies with how v and w are to be read (in some cases they need to be read ulta)
# and then how they need to be printed (alignment_1 vs alignment 2 or ulti)
# It's a hit or miss. I know this works but just not within the file

# NOTE: YOU MAY HAVE TO RUN THIS FROM THE BI3 FOLDER WITH THE COMMAND:
# python -m Wk3.Wk3_3_Linear_Space_Alignment (NO .py)
#########################################################################################
#########################################################################################
# Pseudocode
# Should look to use Wk3_2_MiddleEdge_Linear file here for the midEdge below
"""
LinearSpaceAlignment(v, w, top, bottom, left, right)
    if left = right
        output path formed by bottom − top vertical edges
    if top = bottom
        output path formed by right − left horizontal edges
    middle ← ⌊ (left + right)/2⌋
    midEdge ← MiddleEdge(v, w, top, bottom, left, right)
    midNode ← vertical coordinate of the initial node of midEdge 
    LinearSpaceAlignment(v, w, top, midNode, left, middle)
    output midEdge
    if midEdge = "→" or midEdge = "↘"
        middle ← middle + 1
    if midEdge = "↓" or midEdge ="↘"
        midNode ← midNode + 1
    LinearSpaceAlignment(v, w, midNode, bottom, middle, right)
"""
#########################################################################################


from Wk3.Wk3_2_MiddleEdge_Linear import middle_edge

#########################################################################################
#########################################################################################
# Main Funnction
#########################################################################################
#########################################################################################
def linear_space_alignment(reward, mismatch, indel, v, w):
    """
    Compute the optimal global alignment of two sequences v and w using linear space.

    Args:
        reward: Score for matching characters
        mismatch: Penalty for mismatching characters
        indel: Penalty for insertions/deletions
        v: First sequence (corresponds to rows)
        w: Second sequence (corresponds to columns)
    Returns:
        Tuple containing:
            - Maximum alignment score
            - Aligned version of sequence v
            - Aligned version of sequence w
    """
    # Initialize boundaries for the full matrix
    top, bottom = 0, len(v)
    left, right = 0, len(w)

    # Call the recursive function to compute the alignment
    score, aligned_v, aligned_w = linear_space_alignment_helper(reward, mismatch, indel, v, w, top, bottom, left, right)

    return score, aligned_v, aligned_w

#########################################################################################
#########################################################################################
# Helper functions
#########################################################################################
#########################################################################################
def linear_space_alignment_helper(reward, mismatch, indel, v, w, top, bottom, left, right):
    """
    Helper function for linear space alignment that performs the recursive divide-and-conquer approach.

    Args:
        reward: Score for matching characters
        mismatch: Penalty for mismatching characters
        indel: Penalty for insertions/deletions
        v: First sequence (corresponds to rows)
        w: Second sequence (corresponds to columns)
        top: Top boundary of the current submatrix
        bottom: Bottom boundary of the current submatrix
        left: Left boundary of the current submatrix
        right: Right boundary of the current submatrix
    Returns:
        Tuple containing:
            - Alignment score for the current submatrix
            - Aligned version of the subsequence of v corresponding to the current submatrix
            - Aligned version of the subsequence of w corresponding to the current submatrix
    """
    # Base case: If we are down to a single column, align vertically
    if left == right:
        aligned_v = v[top:bottom]
        aligned_w = '-' * (bottom - top)
        score = -(bottom - top) * indel
        return score, aligned_v, aligned_w

    # Base case: If we are down to a single row, align horizontally
    if top == bottom:
        aligned_v = '-' * (right - left)
        aligned_w = w[left:right]
        score = -(right - left) * indel
        return score, aligned_v, aligned_w

    # Recursive case: Find the middle edge and split the problem
    # middle_edge returns two tuples: (row_start, col_start) and (row_end, col_end)
    # The coordinates are relative to the current submatrix
    mid_edge_start, mid_edge_end = middle_edge(reward, mismatch, indel, v[top:bottom], w[left:right])
    # print('-'*30)
    # print(f"{mid_edge_start=}")
    # print(f"{mid_edge_end=}")


    # Calculate the coordinates of the middle edge in the original matrix
    start_row = top + mid_edge_start[0]
    start_col = left + mid_edge_start[1]
    end_row = top + mid_edge_end[0]
    end_col = left + mid_edge_end[1]

    # print(f"{start_row=}")
    # print(f"{start_col=}")
    # print(f"{end_row=}")
    # print(f"{end_col=}")
    # print('-'*30)

    # Recursively solve for the left half
    left_score, left_aligned_v, left_aligned_w = linear_space_alignment_helper(
        reward, mismatch, indel, v, w, top, start_row, left, start_col)

    # Build the middle edge alignment and calculate its score
    if start_row == end_row:  # Horizontal edge (insertion in v)
        middle_v = '-' * (end_col - start_col)
        middle_w = w[start_col:end_col]
        middle_score = -(end_col - start_col) * indel
    elif start_col == end_col:  # Vertical edge (deletion from v)
        middle_v = v[start_row:end_row]
        middle_w = '-' * (end_row - start_row)
        middle_score = -(end_row - start_row) * indel
    else:  # Diagonal edge
        middle_v = v[start_row]
        middle_w = w[start_col]
        middle_score = reward if middle_v == middle_w else -mismatch

    # Recursively solve for the right half
    right_score, right_aligned_v, right_aligned_w = linear_space_alignment_helper(
        reward, mismatch, indel, v, w, end_row, bottom, end_col, right)

    # Combine all three parts
    aligned_v = left_aligned_v + middle_v + right_aligned_v
    aligned_w = left_aligned_w + middle_w + right_aligned_w
    total_score = left_score + middle_score + right_score

    return total_score, aligned_v, aligned_w



###########################################################################

if __name__ == "__main__":
    # # Sample test
    # reward = 1
    # mismatch = 1
    # indel = 2
    # v = "GAT"
    # w = "GAGA"
    # # Expected answer =
    # # Score: -1
    # # GA-T
    # # GAGA
    # alignment_score, alignment_1, alignment_2 = linear_space_alignment(reward, mismatch, indel, v, w)
    # print(alignment_score)
    # print(alignment_2)
    # print(alignment_1)

    # From file

    # Get dataset
    from pathlib import Path as partho
    current_dir = partho(__file__).parent
    filename = input("Please enter the filename: ")
    file_path = current_dir / filename

    #NOTE: For Wk3_2, there is a newline character after the peptide string, so remove that
    with open(file_path, 'r') as file:
        # Read in n and m
        reward, mismatch, indel = file.readline().split()
        w = file.readline().strip()
        v = file.readline().strip()
        # Down = [list(map(int, line.split())) for line in cleaned_input[1:separator]]

    alignment_score, alignment_1, alignment_2 = linear_space_alignment(int(reward), int(mismatch), int(indel), v, w)
    # print('\n'.join(answer))
    # print(answer)

    with open("Wk3/Wk3_3_output.txt", "w") as output_file:
        output_file.write(str(alignment_score))
        output_file.write("\n")
        output_file.write(alignment_2)
        output_file.write("\n")
        output_file.write(alignment_1)
        # output_file.write(" ".join(map(str, middle_edge_coord_1)))
        # output_file.write("\n")

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