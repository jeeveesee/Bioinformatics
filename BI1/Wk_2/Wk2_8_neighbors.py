##########################################################################################
# TASK:
# Find neighbors of a pattern with up to d mismatches
##########################################################################################


## TESTING:
# text = 'CAGCCTCAGCCTGCAGGAAGCCTCCTGCAGCCTCCTCCTGCAGGCAGCCTCAGCCTCCTGCAGGCAGGAAGAGTAAGTAGCAGGCAGGAAGGCAGCAGGCAGGCAGAGTAAGTAGAAGGAAGAGTAGCAGCAGCAGCCTGAAGAGTACCTCAGCCTGCAGGAAGGCAGAGTAGCAGCAGGCAGCAGGAAGCCTGAAGCAGGCAGCCTCAGCAGGAAGCAGGCAGAGTACCTGAAGCAGCAGCCTCAGGAAGCAGCAGGCAGCCTGAAGGCAGGCAGCAGAGTACAGCCTCAGGCAGGCAGAGTAGAAGAGTACAGCAGGAAGCAGCCTAGTAGAAGAGTAGAAGCCTAGTAAGTAGCAGCCTCAG' 
# k = 5
# d = 3


###############################
# STEP 1: HAMMING DISTANCE
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
# STEP 2: FIND NEIGHBORS
###############################

def neighbors(pattern, d):
    if d == 0:
        return {pattern}
    if len(pattern) == 1:
        return {'A', 'T', 'C', 'G'}

    neighborhood = set()
    suffix_neighbors = neighbors(pattern[1:], d) # Recursive call to get all neighbors without first nucleotide

    for genome in suffix_neighbors:
        if hamming_distance(pattern[1:], genome) < int(d): # Less than d mismatches
            for base in 'ATCG':
                neighborhood.add(base + genome) # Now you get all the neighbors with first base changed
                # The above is the second part of Section 1.7.2 in the textbook
        else:
            neighborhood.add(pattern[0] + genome) 
            # The above keeps the first base remains the same if d mismatches already reached
            # and adds all the genome neighbors to it. This is becase the first base doesn't change,
            # so the hamming distance is still the same.

    return neighborhood


if __name__ == "__main__":
    # Example usage:
    text = input("Please enter your DNA text: ")
    d = input("Please enter the hamming distance, d: ")

    result = neighbors(text, d)
    print("The neighbors are: ", ' '.join(text for text in result))
