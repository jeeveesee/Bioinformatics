#########################################################################################
# TASK:
# For a given list of t DNA strings and an integer k and
# create the best motifs matrix using greedy motif search  with pseudocounts
# Cogniterra Bioinformatics Course: Course 1: Week 3: 1.7 Greedy Motif Search
# GreedyMotifSearch(Dna, k, t)
#     BestMotifs ← motif matrix formed by first k-mers in each string from Dna
#     for each k-mer Motif in the first string from Dna
#         Motif1 ← Motif
#         for i = 2 to t
#             form Profile from motifs Motif1, …, Motifi - 1 using Laplace's rule of succession
#             Motifi ← Profile-most probable k-mer in the i-th string in Dna
#         Motifs ← (Motif1, …, Motift)
#         if Score(Motifs) < Score(BestMotifs)
#             BestMotifs ← Motifs
#     return BestMotifs
# NOTE: NOTE NOTE: NOTE:
# This works but it is not the most efficient. Look at the EFFICIENT FILE for efficiency
##########################################################################################

import numpy as np
import pandas as pd

###############################################
# STEP 1: SET UP MOST PROBABLE K-MER FUNCTION
###############################################

def most_probable_kmer(dna, k, profile_matrix: pd.DataFrame):
    max_prob = -1
    most_prob_kmer = ""

    for i in range(len(dna) - k + 1):
        kmer = dna[i:i+k]
        prob = 1
        # print("kmer: ", kmer)
        # print("Profile_Matrix: \n", profile_matrix)
        for j, nucleotide in enumerate(kmer):
            prob *= profile_matrix.loc[nucleotide, j]
            print(f"nucleotide is {nucleotide} and j is {j}")
            print("Prob so far is: ", prob)

        if prob > max_prob:
            max_prob = prob
            most_prob_kmer = kmer
        
    return max_prob, most_prob_kmer

###############################################
# STEP 2:  GREEDY MOTIF SEARCH w PSEUDOCOUNTS
###############################################


# Function: Create Profile Matrix
# This function creates a profile matrix from a list of motifs
# 1. Start by creating a DF of 1s (ACGT rows, k columns) <--- 1 so that we can add pseudocounts
# 2. Count nucletoide occurrences in each column
# 3. Convert counts to probabilities by dividing by number of motifs or len(motifs list) + 4 (fpor pseudocounts)
#################
def create_profile_matrix(motifs):
    k = len(motifs[0])
    profile = pd.DataFrame(1, index=['A', 'C', 'G', 'T'], columns=range(k), dtype=float)

    # Count occurrences of each nucleotide in each column
    for motif in motifs:
        for j, nucleotide in enumerate(motif):
            profile.at[nucleotide, j] += 1

    # Convert counts to probabilities
    profile = profile.div(len(motifs) + 4)  # +4 for pseudocounts

    return profile


# Function: Score Motifs
# This function creates score of all the nucleotides from all the motifs
# that are not part of the maximum occuring nucleotide in each column
# 1. Start by findding the first element from each motif
# 2. Find the most common nucleotide in this column (this is the most frequent nucleotide in that position)
# 3. Count the remaining nucleotides in that column that are not the most frquent
# 4 Repeat. I am sure that this is done so that we can find the min score 
#################
def score_motifs(motifs):
    k = len(motifs[0])
    t = len(motifs)
    #print(k)
    score = 0
    
    for j in range(k):
        column = [motif[j] for motif in motifs]
        # print(column)
        most_common_nucleotide = max(set(column), key=column.count) 
        # most_frequent_count = max(column.count(nuc) for nuc in 'ATCG')
        score += sum(1 for nucleotide in column if nucleotide != most_common_nucleotide)
        # score += t - most_frequent_count

    return score


def greedy_motif_search_w_pseudocounts(dna_list, k, t):
    # Initialize BestMotifs with the first k-mers from each string in Dna
    best_motifs = [dna[:k] for dna in dna_list]
    # print("Best Motifs: ", best_motifs)
    best_score = score_motifs(best_motifs)
    # print("Best Score: ", best_score)

    first_dna = dna_list[0]
    # print("First DNA: ", first_dna)

    for i in range(len(first_dna) - k + 1):
        motif1 = first_dna[i:i+k]
        # print("Motif 1: ", motif1)
        motifs = [motif1]
        # print("Motifs: ", motifs)

        for j in range(1, t):
            current_profile = create_profile_matrix(motifs)
            # print("Current Profile: \n", current_profile)
            _, next_motif = most_probable_kmer(dna_list[j], k, current_profile)
            # print(f"Next Motif from DNA string {j+1}: ", next_motif)
            motifs.append(next_motif)
            # print("Updated Motifs: ", motifs)

        current_score = score_motifs(motifs)
        # print("Current Score: ", current_score)

        if current_score < best_score:
            best_score = current_score
            best_motifs = motifs

    return best_motifs


# Example from text
print(greedy_motif_search_w_pseudocounts(['GGCGTTCAGGCA', 'AAGAATCAGTCA', 'CAAGGAGTTCGC', 'CACGTCAATCAC', 'CAATAATATTCG'], 3, 5))
# print(create_profile_matrix(['TCATGAGTAGTC']))
# print(create_profile_matrix(['TCA', 'TGA', 'GTA', 'GTC']))
# print(score_motifs(['GGC', 'AAG', 'CAA', 'CAC', 'CAA']))
# print(create_profile_matrix(['TCA', 'TGA', 'GTA', 'GTC']))
# print(score_motifs(['GGC', 'AAG', 'CAA', 'CAC', 'CAA']))
# print(most_probable_kmer('AAGGAGTTCGC', 3, create_profile_matrix(['TCA', 'TGA', 'GTA', 'GTC'])))


# My dataset
# if __name__ == "__main__":
#     filename = input("Please enter the filename: ")

#     with open(filename, 'r') as file:
#             # Read k and t (Row 0)
#             k_and_t = file.readline().strip()
#             k, t = map(int, k_and_t.split())

#             # Read the DNA strings
#             dna_list_raw = file.readline().strip()
#             dna_list = dna_list_raw.split()

#     answer = greedy_motif_search_w_pseudocounts(dna_list, k, t)
#     print(*answer)

