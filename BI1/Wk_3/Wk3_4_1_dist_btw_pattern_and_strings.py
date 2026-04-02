#########################################################################################
# TASK:
# Find 
# d(Pattern, DNA)
# where d(Pattern, DNA) is the sum of the distances between Pattern and each string in DNA
# This seems like a subroutine of the broader Median String Problem given later
# in the text - ah well, let's do it anyway
##########################################################################################

# from itertools import product

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
# STEP 3: RUN THE PROGRAM
###############################


if __name__ == "__main__":
    # Example usage:
    dna_list_raw = input("Please enter your DNA strings separated by a space: ")
    dna_list = dna_list_raw.split()
    #print(dna_list)
    kmer = input("Please enter your k-mer pattern: ")

    dist = pattern_distance(kmer, dna_list)
    print("The distance is:", dist)