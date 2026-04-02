#########################################################################################
# TASK:
# Use brute force method to a k-mer Pattern (or Median String) that minimizes
# d(Pattern, DNA), among all possible choices of k-mers
# where d(Pattern, DNA) is the sum of the distances between Pattern and each string in DNA
##########################################################################################

from itertools import product

###############################
# STEP 1: Hamming Distance
###############################

def hamming_distance(p_genome, q_genome):

    if len(p_genome) != len(q_genome):
        return "Error: Strings must be of equal length"

    hamming_dist = 0
    for i in range(len(p_genome)):
        if p_genome[i] != q_genome[i]:
            hamming_dist += 1

    return hamming_dist


###############################
# STEP 2: Pattern Distance
###############################

def pattern_distance(pattern, dna):
    k = len(pattern)
    total_distance = 0

    for string in dna:
        min_distance = float('inf')
        for i in range(len(string) - k + 1):
            kmer = string[i: i+k]
            current_distance = hamming_distance(pattern, kmer)
            if current_distance < min_distance:
                min_distance = current_distance
        total_distance += min_distance

    return total_distance


###############################
# STEP 3: MEDIAN STRING
###############################

def median_string(dna, k):
    distance = float('inf')
    # Generate all kmer possibilities
    all_kmers = [''.join(v) for v in product('ATCG', repeat=k)]
    
    for kmer in all_kmers:
        current_distance = pattern_distance(kmer, dna)
        if current_distance < distance:
            distance = current_distance
            median = kmer

    return median


###############################
# STEP 3: RUN THE PROGRAM
###############################


if __name__ == "__main__":
    # Example usage:
    dna_list_raw = input("Please enter your DNA strings separated by a space: ")
    dna_list = dna_list_raw.split()
    print(dna_list)
    k = input("Please enter the k-mer length, k: ")

    median = median_string(dna_list, int(k))
    print("The median string is:", median)