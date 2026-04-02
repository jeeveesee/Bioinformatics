#########################################################################################################################################
# TASK:
# Use Gibbs Sampler method to find motifs that might be the
# binding site for the transcription factor regulating these genes 

"""
Mycobacterium tuberculosis (MTB) can persist in a latent state in humans 
for many years before causing disease. 
Latency has been found to be linked to hypoxia (lack  of oxygen) in the host. 
You suspect that genes that are activated in  hypoxia are regulated by a 
common transcription factor, so you collect  the upstream sequences for all 
of the MTB genes that are upregulated in  hypoxia, 
looking for the motif that corresponds to the binding site for  the 
transcription factor regulating these genes. Your biologist  
colleague tells you that you should look at the 250 bp upstream region  of 
each gene (which have been conveniently compiled for you in a FASTA
file named upstream250.txt right click and download this file). 
Your colleague also tells you that the motif is probably about 20 bp long.
"""
#########################################################################################################################################

# Imports
import random

# Main Gibbs Function
#########################################################################################################################################

# Gibbs Sampler Function
def gibbs_sampler(dna_list, k, t, N, num_restarts=20):
    """ Utilizes Gibbs Sampling and randomly generated kmers with resubstitution to find the best concensus motifs in a DNA list"""

    best_overall_motifs = None
    best_overall_score = float('inf')

    for _ in range(num_restarts):

        #Initialize random kmer motifs from each dna string
        motifs = []
        for dna in dna_list:
            start_index = random.randint(0, len(dna) - k)
            motifs.append(dna[start_index:start_index + k])

        best_motifs = motifs[:]

        # Run Gibbs Sampler with resubstitution N times
        for _ in range(N):
            i = random.randint(0, t-1)
            motifs_except_i = motifs[:i] + motifs[i+1:]
            profile = create_profile_matrix(motifs_except_i)
            new_motif = most_probable_kmer(dna_list[i], k, profile)
            motifs[i] = new_motif
            if score_motifs(motifs) < score_motifs(best_motifs):
                best_motifs = motifs[:]

        current_score = score_motifs(best_motifs)
        if current_score < best_overall_score:
            best_overall_score = current_score
            best_overall_motifs = best_motifs[:]

    return best_overall_motifs


# Helper Functions
#########################################################################################################################################

# 1: Profile Matrix
#########################
def create_profile_matrix(motifs):
    """ Creates a profile matrix from a list of motifs with pseudocounts"""
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
        most_frequent_count = max(column.count(nuc) for nuc in 'ACGT')
        score += t - most_frequent_count # How many nucleotide are NOT the most frequent in this column

    return score

# 4. Form consensus string from motifs
###########################################################

def form_consensus_string(motifs):
    """ Forms a consensus string from a list of motifs"""
    k = len(motifs[0])
    consensus = []

    for j in range(k):
        column = [motif[j] for motif in motifs]
        most_frequent_nuc = max('ACGT', key=lambda nuc: column.count(nuc))
        consensus.append(most_frequent_nuc)

    return ''.join(consensus)

# Run Example / Exam
#########################################################################################################################################

if __name__ == "__main__":

    # Example from text
    # dna_list = ['CGCCCCTCTCGGGGGTGTTCAGTAACCGGCCA', 'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG', 'TAGTACCGAGACCGAAAGAAGTATACAGGCGT', 'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC', 'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']
    # k = 8
    # t = 5
    # N = 100

    # best_answer = gibbs_sampler(dna_list, k, t, N, num_restarts=20)
    # print("\n\nThe best motifs are: ", *best_answer)

    # Answer: TCTCGGGG CCAAGGTG TACAGGCG TTCAGGTG TCCACGTG

   #random.seed(42)  # For reproducibility <--- THIS IS WHAT WAS MESSING EVERYTHING UP!! DON'T USE A SEED
   # Get dataset

# Uncomment the below for Gibbs Sampler Test Set
   from pathlib import Path as partho
   current_dir = partho(__file__).parent
   filename = input("Please enter the filename: ")
   file_path = current_dir / filename

   dna_list = []
   with open(file_path, 'r') as file:

       #Read sequences except for those starting with ">"
       for line in file:
           if not line.startswith('>'):
               dna_list.append(line.strip())

       num_restarts = 20
       k = 20
       t = len(dna_list)
       N = 2000

       best_answer = gibbs_sampler(dna_list, k, t, N, num_restarts=num_restarts)
       print("\n\nThe best motifs are: ", *best_answer)
       print("\nThe consensus string is: ", form_consensus_string(best_answer))



