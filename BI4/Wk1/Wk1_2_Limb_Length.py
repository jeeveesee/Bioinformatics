##############################################################################################
# Molecular Evolution  - Bioinformatic IV Course from Coursera
#
# Week 1 - Distance-Based Phylogeny Construction
# Code Challenge: Solve the Limb Length Problem.
# Input: An integer n, followed by an integer j between 0 and n - 1,
# followed by a space-separated additive distance matrix D (whose elements are integers).
# Output: The limb length of the leaf in Tree(D) corresponding to row j of this distance matrix
# (use 0-based indexing).
##############################################################################################
"""
Pseudocode (Compeau & Pevzner):

LimbLength(j, n, D)
    for each pair of leaves i and k such that i ≠ j and k ≠ j
        limb_length(j) <- (D(i, j) + D(j, k) - D(i, k)) / 2
    return minimum limb_length(j) over all valid pairs i, k
"""

from pathlib import Path as partho


# Main code

def parse_distance_matrix(dist_matrix_str):
    """
    Parse a tab/space-separated distance matrix string into a 2D list of integers.

    Parameters:
        dist_matrix_str (str): Multi-line string with tab or space-separated integer values.

    Returns:
        list[list[int]]: 2D list representing the distance matrix.
    """
    matrix = []
    for line in dist_matrix_str.strip().split('\n'):
        row = list(map(int, line.split()))
        matrix.append(row)
    return matrix


def limb_length(n, j, dist_matrix_str):
    """
    Compute the limb length of leaf j in the tree represented by an additive distance matrix.

    Uses the formula:
        LimbLength(j) = min over all pairs i, k (i != j, k != j, i != k)
                        of (D[i][j] + D[j][k] - D[i][k]) / 2

    Parameters:
        n (int): Number of leaves in the tree (used for input parsing context).
        j (int): 0-based index of the leaf whose limb length is to be computed.
        dist_matrix_str (str): Multi-line string of the additive distance matrix.

    Returns:
        int: The limb length of leaf j.
    """
    d = parse_distance_matrix(dist_matrix_str)
    # print(f"{d=}")
    size = len(d) # Num of rows
    # print(f"{size=}")

    min_limb = float('inf')
    for i in range(size):
        # print(f"Orig {i=}")
        # print(f"Orig {j=}")
        if i == j:
            continue
        for k in range(i + 1, size):
            # print(f"{i=}")
            # print(f"{j=}")
            # print(f"{k=}")
            # print(f"{min_limb=}")
            if k == j:
                continue
            limb = (d[i][j] + d[j][k] - d[i][k]) // 2
            # print(f"  i={i}, j={j}, k={k}: limb = ({d[i][j]} + {d[j][k]} - {d[i][k]}) // 2 = {limb}")
            if limb < min_limb:
                min_limb = limb

    return min_limb


def formatterer(answer):
    """
    Format the limb length result as a plain integer string.

    Parameters:
        answer (int): The computed limb length.

    Returns:
        str: The answer as a string.
    """
    return str(answer)


###########################################################################

if __name__ == "__main__":
    # Sample test
    n = 4
    j = 1
    dist_matrix_str = """0	13	21	22
13	0	12	13
21	12	0	13
22	13	13	0"""
    # Expected answer =
    # 2
    answer = limb_length(n, j, dist_matrix_str)
    print(formatterer(answer)) # Formatted answer

    # # From file

    # # Get dataset
    # from pathlib import Path as partho

    # current_dir = partho(__file__).parent
    # filename = input("Please enter the filename: ")
    # file_path = current_dir / filename

    # with open(file_path, "r") as file:
    #     data = file.read().strip()
    #     lines = data.split('\n')
    #     n = int(lines[0])
    #     j = int(lines[1])
    #     dist_matrix_str = '\n'.join(lines[2:])

    # answer = limb_length(n, j, dist_matrix_str)

    # with open("Wk1_2_output.txt", "w") as output_file:
    #     output_file.write(formatterer(answer))
