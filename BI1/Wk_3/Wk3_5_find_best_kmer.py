#########################################################################################
# TASK:
# For a given text string, an integerk and a 4xk matrix profil
# find the most probably k-mer in the string
# Basically, use greedy approach to find a k-mer, then multiply the probabilities
##########################################################################################

import numpy as np
import pandas as pd

def most_probable_kmer(dna, k, profile_matrix):
    max_prob = -1
    most_prob_kmer = ""

    for i in range(len(dna) - k + 1):
        kmer = dna[i:i+k]
        # print("kmer: ", kmer)
        # print("Profile_Matrix: \n", profile_matrix)
        prob = 1
        
        for j, nucleotide in enumerate(kmer):
            # print(nucleotide, j)
            prob *= profile_matrix.loc[nucleotide, j]
            # print(f"nucleotide is {nucleotide} and j is {j}")
            # print("Prob so far is: ", prob)
            ## Crappy Way
            ## if nucleotide == 'A':
            ##     prob *= profile_matrix[0][j]
            ##     #print("A", profile_matrix[0][j])
            ## elif nucleotide == 'C':
            ##     prob *= profile_matrix[1][j]
            ##     #print("C", profile_matrix[1][j])
            ## elif nucleotide == 'G':
            ##     prob *= profile_matrix[2][j]
            ## elif nucleotide == 'T':
            ##     prob *= profile_matrix[3][j]
            ##print(nucleotide, j, prob, max_prob)
        
        if prob > max_prob:
            max_prob = prob
            most_prob_kmer = kmer
        
    return max_prob, most_prob_kmer


if __name__ == "__main__":
    # Example usage:
    # dna_raw = input("Please enter your DNA string: ")
    # dna = dna_raw.split()
    # k = input("Please enter the k-mer length, k: ")
    # profile_matrix = np.array([
    #                     [0.2, 0.2, 0.3, 0.2, 0.3],
    #                     [0.4, 0.3, 0.1, 0.5, 0.1],
    #                     [0.3, 0.3, 0.5, 0.2, 0.4],
    #                     [0.1, 0.2, 0.1, 0.1, 0.2]
    # ])
    filename = input("Please enter the filename: ")

    with open(filename, 'r') as file:
        # Read the DNA string (Row 0)
        dna = file.readline().strip()

        # Read k (Row 1) and convert to integer
        k_str = file.readline().strip()
        try:
            k = int(k_str)
        except ValueError:
            raise ValueError(f"Could not convert '{k_str}' to an integer for k.")

        # Read the remaining lines for the profile matrix (Rows 2-5)
        profile_matrix_np = np.genfromtxt(file, dtype=float)
        profile_matrix = pd.DataFrame(profile_matrix_np, index = ['A', 'C', 'G', 'T'])

    max_prob_ans, max_kmer_ans = most_probable_kmer(dna, int(k), profile_matrix)
    print(f"Max Probability is: {max_prob_ans} and the max kmer is: {max_kmer_ans}" )