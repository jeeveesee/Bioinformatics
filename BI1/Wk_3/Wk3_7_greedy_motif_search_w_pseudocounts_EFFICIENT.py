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
##########################################################################################

###############################################
# STEP 1: SET UP MOST PROBABLE K-MER FUNCTION
###############################################

def most_probable_kmer(dna, k, profile_matrix):
    max_prob = -1
    most_prob_kmer = ""

    for i in range(len(dna) - k + 1):
        kmer = dna[i:i+k]
        prob = 1
        for j, nucleotide in enumerate(kmer):
            prob *= profile_matrix[nucleotide][j]

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
# 3. Convert counts to probabilities by dividing by number of motifs or len(motifs list) + 4 (for pseudocounts)
# Doing this without pandas for efficiency
#################
def create_profile_matrix(motifs):
    k = len(motifs[0])
    t = len(motifs)
    pseudocount = 1
    profile = {nuc: [pseudocount] * k for nuc in 'ACGT'} # Form profile of 1s for pseudocounts for each of ATGC

    # Count occurrences of each nucleotide in each column
    for motif in motifs:
        for j, nucleotide in enumerate(motif):
            profile[nucleotide][j] += 1

    # Convert counts to probabilities
    # +4 for pseudocounts since we added 1 to each of the 4 nucleotides
    profile = {key: [value/(t+4) for value in values]
               for key, values in profile.items()}

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
    score = 0

    for j in range(k):
        column = [motif[j] for motif in motifs]
        # Don't need to find the nucleotide, just its count
        most_frequent_count = max(column.count(nuc) for nuc in 'ATCG')
        score += t - most_frequent_count # How many nucleotide are NOT the most frequent in this column

    return score

# Final greedy search
#################
def greedy_motif_search_w_pseudocounts(dna_list, k, t):
    # Initialize BestMotifs with the first k-mers from each string in Dna
    best_motifs = [dna[:k] for dna in dna_list]
    best_score = score_motifs(best_motifs)
    # print("Best Score: ", best_score)
    first_dna = dna_list[0]

    for i in range(len(first_dna) - k + 1):
        motif1 = first_dna[i:i+k]
        motifs = [motif1]

        for j in range(1, t):
            current_profile = create_profile_matrix(motifs)
            _, next_motif = most_probable_kmer(dna_list[j], k, current_profile)
            motifs.append(next_motif)

        current_score = score_motifs(motifs)

        if current_score < best_score:
            best_score = current_score
            best_motifs = motifs

    return best_motifs


#################
# Example from text
# print(greedy_motif_search_w_pseudocounts(['GGCGTTCAGGCA', 'AAGAATCAGTCA', 'CAAGGAGTTCGC', 'CACGTCAATCAC', 'CAATAATATTCG'], 3, 5))
# ANSWER: ['TTC', 'ATC', 'TTC', 'ATC', 'TTC']
# print(create_profile_matrix(['TCATGAGTAGTC']))
# print(create_profile_matrix(['TCA', 'TGA', 'GTA', 'GTC']))
# print(score_motifs(['GGC', 'AAG', 'CAA', 'CAC', 'CAA']))
# print(most_probable_kmer('AAGGAGTTCGC', 3, {'A': [0.125, 0.125, 0.5], 'C': [0.125, 0.25, 0.25], 'G': [0.375, 0.25, 0.125], 'T': [0.375, 0.375, 0.125]}))

# My dataset
if __name__ == "__main__":
    filename = input("Please enter the filename: ")

    with open(filename, 'r') as file:
            # Read k and t (Row 0)
            k_and_t = file.readline().strip()
            k, t = map(int, k_and_t.split())

            # Read the DNA strings
            dna_list_raw = file.readline().strip()
            dna_list = dna_list_raw.split()

    answer = greedy_motif_search_w_pseudocounts(dna_list, k, t)
    print(*answer)

