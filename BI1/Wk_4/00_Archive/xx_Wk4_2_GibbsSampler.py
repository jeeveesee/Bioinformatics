#########################################################################################################################################
# TASK:
# Use Gibbs Sampler method to generate  profiles to
# find a set of motifs that are likely implanted by a common motif
# NOTE: YOU HAVE TO RUN THIS FROM THE BI1 FOLDER WITH THE COMMAND:
# python -m Wk_4.Wk4_2_GibbsSampler

"""
GibbsSampler(Dna, k, t, N)
    randomly select k-mers Motifs = (Motif1, …, Motift) in each string from Dna
    BestMotifs ← Motifs
    for j ← 1 to N
        i ← Random(t)
        Profile ← profile matrix constructed from all strings in Motifs except for Motifi
        Motifi ← Profile-randomly generated k-mer in the i-th sequence
        if Score(Motifs) < Score(BestMotifs)
            BestMotifs ← Motifs
    return BestMotifs
"""
#########################################################################################################################################

# Imports
import random

# Main Gibbs Function
#########################################################################################################################################

# Gibbs Sampler Function
def gibbs_sampler(dna_list, k, t, n):
    """ Utilizes Gibbs Sampling and randomly generated kmers with resubstitution to find the best concensus motifs in a DNA list"""
    # Initialize random kmer motifs from each DNA string
    motifs = []
    for dna in dna_list:
        start_index = random.randint(0, len(dna) - k)
        motifs.append(dna[start_index:start_index + k])

    best_motifs = motifs[:]

    # Run Gibbs Sampler with resubstitution N times
    for _ in range(n):
        i = random.randint(0, t-1)
        motifs_except_i = motifs[:i] + motifs[i+1:]
        profile = create_profile_matrix(motifs_except_i)
        new_motif = most_probable_kmer(dna_list[i], k, profile)
        motifs[i] = new_motif
        if score_motifs(motifs) < score_motifs(best_motifs):
            best_motifs = motifs[:]

    return best_motifs

# Helper Functions
#########################################################################################################################################

# 1: Profile Matrix
#########################
def create_profile_matrix(motifs):
    """ Creates a profile matrix from a list of motifs with pseudocounts"""
    k = len(motifs[0])
    t = len(motifs)
    pseudocount = 1
    profile = {nuc: [pseudocount] * k for nuc in 'ATCG'} # Form profile of 1s for pseudocounts for each of ATGC

    # Count occurrences of each nucleotide in each column
    for motif in motifs:
        for j, nucleotide in enumerate(motif):
            profile[nucleotide][j] += 1

    # Convert counts to probabilities
    # +4 for pseudocounts since we added 1 to each of the 4 nucleotides
    profile = {key: [value/(t+4) for value in values]
               for key, values in profile.items()}

    return profile

# 2. Most Probable kmer Function with Pseudocounts and biased profile
########################################################################
def most_probable_kmer(dna, k, profile_matrix):
    """ Provides most problable kmers in a DNA string giveb a profile matrix with pseudocounts"""

    kmers = [] # To store all the kmers
    probs = [] # To store all probs that will then be used to randomly select a kmer with weighted probaility

    for i in range(len(dna) - k + 1):
        kmer = dna[i:i+k]
        prob = 1
        for j, nucleotide in enumerate(kmer):
            prob *= profile_matrix[nucleotide][j]

        probs.append(prob)
        kmers.append(kmer)

    # Normalize probablities to sum to 1
    total_prob = sum(probs)
    normalized_probs = [p / total_prob for p in probs]

    # Randomly select kmer based on weighted probabilities
    most_prob_kmer = random.choices(kmers, weights=normalized_probs, k=1)[0]

    return most_prob_kmer


# 3. Score Motifs
###########################################################
def score_motifs(motifs):
    """ Scores motifs by summing the number of nucleotides not equal to the most freq nucleotide in each column"""
    k = len(motifs[0])
    t = len(motifs)
    score = 0

    for j in range(k):
        column = [motif[j] for motif in motifs]
        # Don't need to find the nucleotide, just its count
        most_frequent_count = max(column.count(nuc) for nuc in 'ATCG')
        score += t - most_frequent_count # How many nucleotide are NOT the most frequent in this column

    return score


# Run Example / Exam
#########################################################################################################################################

if __name__ == "__main__":
    
    # Example from text
    dna_list = ['CGCCCCTCTCGGGGGTGTTCAGTAACCGGCCA', 'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG', 'TAGTACCGAGACCGAAAGAAGTATACAGGCGT', 'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC', 'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']
    k = 8
    t = 5
    n = 5000

    best_answer = gibbs_sampler(dna_list, k, t, n)
    print("\n\nThe best motifs are: ", *best_answer)

    # Answer: TCTCGGGG CCAAGGTG TACAGGCG TTCAGGTG TCCACGTG

   #random.seed(42)  # For reproducibility <--- THIS IS WHAT WAS MESSING EVERYTHING UP!! DON'T USE A SEED
   # Get dataset

# Uncomment the below for Gibbs Sampler Test Set
#    from pathlib import Path as partho
#    current_dir = partho(__file__).parent
#    filename = input("Please enter the filename: ")
#    file_path = current_dir / filename

#    with open(file_path, 'r') as file:

#        # Read k, t and N (Row 0)
#        k_t_N = file.readline().strip()
#        k, t, N = map(int, k_t_N.split())

#        # Read the DNA strings
#        dna_list_raw = file.readline().strip()
#        dna_list = dna_list_raw.split()

    #    best_answer = gibbs_sampler(dna_list, k, t, N)
    #    print("\n\nThe best motifs are: ", *best_answer)



