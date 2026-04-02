##########################################################################################
# TASK:
# Find approximate match of a kmer in a string
# Utilize the hamming distance to find an instande where the kmer
# matches within a certain distance
# Add that to a list
##########################################################################################

###############################
# HAMMING DISTANCE
###############################

def hamming_distance(p_genome, q_genome):
    
    if len(p_genome) != len(q_genome):
        return "Error: Strings must be of equal length"
    
    hamming_dist = 0
    for i in range(len(p_genome)):
        if p_genome[i] != q_genome[i]:
            hamming_dist += 1
    return hamming_dist


#################################
# APPROXIMATE MATCHING PROBLEM
#################################

# Testing
# text = "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT"
# kmer = "ATTCTGGA"
# dist = 3
# text = "TTTAGAGCCTTCAGAGG"
# kmer = "GAGG"
# dist = 2

# READ FILE
file_name = input("Please enter file name: ")
with open(file_name, 'r') as file:
    content = file.read()
    genomes = content.splitlines()
    kmer = genomes[0]
    text = genomes[1]
    dist = genomes[2]

# # Quiz
# text = 'TACGCATTACAAAGCACA'
# kmer = 'AA'
# dist = 1

def approx_match(text, kmer, dist):
    approx_match_list = []
    for i in range(len(text) - len(kmer)+1):
        p_genome = kmer
        q_genome = text[i:i+len(kmer)]

        if(hamming_distance(p_genome, q_genome)) <= int(dist):
            approx_match_list.append(i)
    
    return approx_match_list

print("The approximate matches are at indices: ", ' '.join(str(item) for item in approx_match(text, kmer, dist)))
print("\n\n")
print(f"Total Matches: {len(approx_match(text, kmer, dist))}")
