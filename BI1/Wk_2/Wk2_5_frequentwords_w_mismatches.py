##########################################################################################
# TASK:
# Find frequent words of length k with up to d mismatches in Text
# Output: all most frequent k-mers with up to d mismatches in Text
# Create a pattern dictionary with all possible k-mers and their neighbors
# Find the maximum value in the dictionary
# Add that to a list

"""
Text from textbook:
For example, to generate Neighbors(CAA,1), 
first form Neighbors(AA,1) = {AA, CA, GA, TA, AC, AG, AT}. 
The Hamming distance between AA and each of six of these neighbors is 1. 
Firstly, concatenating C with each of these patterns results in six patterns 
(CAA, CCA, CGA, CTA, CAC, CAG, CAT) that belong to Neighbors(CAA, 1). 
Secondly, concatenating any nucleotide with AA (of CAA) results in four patterns 
(AAA, CAA, GAA, and TAA) that belong to Neighbors(CAA, 1). 
Thus, Neighbors(CAA, 1) comprises ten patterns.
"""

"""
Freq words with mismatches:
One way to solve the Frequent Words with Mismatches problem is to 
generate all 4k k-mers Pattern, 
compute ApproximatePatternCount(Text, Pattern, d) for each k-mer Pattern, 
and then find k-mers with the maximum number of approximate occurrences. 
This is an inefficient approach in practice, 
since many of the 4k k-mers should not be considered 
because neither they nor their mutated versions (with up to d mismatches) 
appear in Text.

Instead, the following pseudocode will generalize the BetterFrequentWords function 
and its use of the frequency table. 
It uses a single map that counts the number of times a given 
string has an approximate match in Text. 
For a given k-mer substring Pattern of Text, 
we need to increase 1 to the count of every k-mer 
that has Hamming distance at most d from Pattern.  
The collection of all such k-mers is called the 
d-neighborhood of Pattern, denoted Neighbors(Pattern, d).
"""
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
        if hamming_distance(pattern[1:], genome) < d: # Less than d mismatches
            for base in 'ATCG':
                neighborhood.add(base + genome) # Now you get all the neighbors with first base changed
                # The above is the second part of Section 1.7.2 in the textbook
        else:
            neighborhood.add(pattern[0] + genome) 
            # The above keeps the first base remains the same if d mismatches already reached
            # and adds all the genome neighbors to it. This is becase the first base doesn't change,
            # so the hamming distance is still the same.

    return neighborhood



########################################
# STEP 3: FREQ WORDS WITH MISMATCHES
########################################

def frequent_words_with_mismatches(text, k, d):
    patterns_dict = {}
    n = len(text)
    
    for i in range(n - k + 1):
        pattern = text[i:i+k]
        neighborhood = neighbors(pattern, d)
    
        for approx_pattern in neighborhood:
            if approx_pattern in patterns_dict:
                patterns_dict[approx_pattern] += 1
            else:
                patterns_dict[approx_pattern] = 1
    
    max_pattern_value = max(patterns_dict.values())
    frequent_patterns = [pattern_key\
                         for pattern_key, value_count in patterns_dict.items()\
                         if value_count == max_pattern_value]

    return frequent_patterns, max_pattern_value

text = input("Please enter your DNA text: ")
k = input("Please enter the length of k-mer, k: ")
d = input("Please enter the hamming distance, d: ")

result = frequent_words_with_mismatches(text, int(k), int(d))
print("The most frequent k-mers with up to d mismatches are: ", ' '.join(result[0]))
print("The maximum number of occurrences is: ", result[1])
print("The number of k-mers with maximum occurrences is: ", len(result[0]))
