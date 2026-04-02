# from Wk_2.Wk2_8_neighbors import hamming_distance, neighbors

# p_genome = 'ATGC'
# q_genome = 'ATCG'

# print(hamming_distance(p_genome, q_genome))

# Read in file

import numpy as np

def read_data_from_file(filepath):
    """
    Reads DNA string, k, and a profile matrix from a file.

    Args:
        filepath (str): The path to the input file.

    Returns:
        tuple: (dna_string, k_integer, profile_matrix_numpy_array)
    """
    with open('data.txt', 'r') as file:
        # Read the DNA string (Row 0)
        dna = file.readline().strip()

        # Read k (Row 1) and convert to integer
        k_str = file.readline().strip()
        try:
            k = int(k_str)
        except ValueError:
            raise ValueError(f"Could not convert '{k_str}' to an integer for k.")


        # Read the remaining lines for the profile matrix (Rows 2-5)
        # and convert them into a list of lists of floats
        # matrix_data = []
        # for line in file:
        #     # Split the line by spaces, filter out empty strings (in case of extra spaces),
        #     # and convert to float
        #     row = [float(x) for x in line.strip().split() if x]
        #     if row: # Only append non-empty rows
        #         matrix_data.append(row)

        profile_matrix = np.genfromtxt(file, dtype=float)

    # # Convert the list of lists into a NumPy array
    # profile_matrix = np.array(matrix_data)

    return dna, k, profile_matrix

print(read_data_from_file('data.txt')[0])
print(read_data_from_file('data.txt')[1])
print(read_data_from_file('data.txt')[2])