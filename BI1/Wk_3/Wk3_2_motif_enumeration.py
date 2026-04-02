#########################################################################################
# TASK:
# Use brute force method to find all k-mers within a d mismatch 
# of a collection of DNA strings
##########################################################################################

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

###############################
# STEP 3: MOTIF ENUMERATION
###############################

def motif_enumeration(dna_list, k, d):

    patterns = set()
    first_dna_string = dna_list[0]
    n = len(first_dna_string)

    for i in range(n-k+1): # Run through all k-mers from the first string
        kmer = first_dna_string[i:i+k]
        neighborhood = neighbors(kmer, d)
        
        for neighbor in neighborhood: # Look for neighbors in all other strings
            found_in_all = True # Assume it's found in all strings unless proven otherwise

            for string in dna_list[1:]: # Search for each neighbor in each string
                found_in_current = False # Assume it's not found in the current string unless proven otherwise

                for j in range(len(string)-k+1): # Searching all k-mers in THIS string
                    current_kmer = string[j:j+k]
                    if hamming_distance(neighbor, current_kmer) <= int(d):
                        found_in_current = True
                        # breaking here means we found neighbor in current string
                        # no need to keep looking in this string, move to next string
                        break 
                
                # If not found in current string, break and move to next neighbor 
                # as this neighbor cannot be in all strings
                if not found_in_current: 
                    found_in_all = False
                    break
            
            # if found in all strings, add to patterns set
            if found_in_all:
                patterns.add(neighbor)
    
    return patterns


if __name__ == "__main__":
    # Example usage:
    dna_list_raw = input("Please enter your DNA strings separated by a space: ")
    dna_list = dna_list_raw.split()
    k = input("Please enter the k-mer length, k: ")
    d = input("Please enter the hamming distance, d: ")

    result = motif_enumeration(dna_list, int(k), int(d))
    print("The motifs are: ", ' '.join(text for text in result))
