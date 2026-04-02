#########################################################################################
# Genome Sequencing  - Bioinformatic III Course from Coursera
#
# Week 2 -Levenshtein Edit distance
# Edit Distance Problem: Find the edit distance between two strings.
# Input: Two strings.
# Output: The edit distance between these strings.
# NOTE: YOU HAVE TO RUN THIS FROM THE BI3 FOLDER WITH THE COMMAND:
# python -m Wk2.Wk2_2_Local_Alignment (NO .py)
#########################################################################################

def edit_distance(v, w):
    """
    Uses Levenshtein Edit Distance to find the minimum number of edit operations required
    to convert string v to string w

    Parameters:
    v, w -> two sequences

    Returns:
    Edit distance between the two strings
    """
    # Initialize scoring matrix
    n = len(v)
    m = len(w)
    S = [[0 for j in range(m + 1)] for i in range(n + 1)]
    # print(f"Pre-initialization, {S=}")

    # Initialize first row and column
    for i in range(n + 1):
        S[i][0] = i
    for j in range(m + 1):
        S[0][j] = j
    # print(f"Post-initialization, {S=}")

    # Fill in the scoring matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # print("")
            # print(f"{i=}")
            # print(f"{i-1=}")
            # print(f"{j=}")
            # print(f"{j-1=}")
            # print(f"{v[i-1]=}")
            # print(f"{w[j-1]=}")
            if v[i - 1] == w[j - 1]:
                # print("Yes, they match")
                cost = 0
            else:
                cost = 1
                # print("No, they don't match")
            # print(f"{S[i-1][j]=}")
            # print(f"{S[i][j-1]=}")
            # print(f"{S[i-1][j-1]=}")
            S[i][j] = min(S[i - 1][j] + 1,      # Deletion
                          S[i][j - 1] + 1,      # Insertion
                          S[i - 1][j - 1] + cost)  # Match/Mismatch
            # print(f"{S[i][j]=}")

    return S[n][m]

###########################################################################

if __name__ == "__main__":

    # # Sample test
    # v = "GAGA"
    # w = "GAT"
    # # Expected answer =
    # # 2
    # answer = edit_distance(v, w)
    # # print('\n'.join(answer))
    # print(answer)

# From file

    # Get dataset
    from pathlib import Path as partho
    current_dir = partho(__file__).parent
    filename = input("Please enter the filename: ")
    file_path = current_dir / filename

    #NOTE: For Wk3_2, there is a newline character after the peptide string, so remove that
    with open(file_path, 'r') as file:
        # reward, mismatch, indel_sigma = file.readline().split()
        v = file.readline().strip()
        w = file.readline().strip()
        # Down = [list(map(int, line.split())) for line in cleaned_input[1:separator]]

    answer = edit_distance(v, w)
    # print('\n'.join(answer))
    print(str(answer))

    with open("Wk2_3_output.txt", "w") as output_file:
        output_file.write(str(answer))


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