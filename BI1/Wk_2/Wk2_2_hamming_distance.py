##########################################################################################
# TASK:
# Create a function to find the hamming distance between two strings
# if p[i] != q[i] then increase the hamming distance by 1
##########################################################################################

# Read file
file_name = input("Please enter file name: ")
with open(file_name, 'r') as file:
    content = file.read()
    genomes = content.splitlines()
    p_genome = genomes[0]
    q_genome = genomes[1]
    #print(f"Genome 1: {p_genome}")
    #print(f"Genome 2: {q_genome}")

# # Quiz
# p_genome = 'CTTGAAGTGGACCTCTAGTTCCTCTACAAAGAACAGGTTGACCTGTCGCGAAG'
# q_genome = 'ATGCCTTACCTAGATGCAATGACGGACGTATTCCTTTTGCCTCAACGGCTCCT'

def hamming_distance(p_genome, q_genome):
    
    if len(p_genome) != len(q_genome):
        return "Error: Strings must be of equal length"
    
    hamming_dist = 0
    for i in range(len(p_genome)):
        if p_genome[i] != q_genome[i]:
            hamming_dist += 1
    return hamming_dist

print(f"Hamming Distance: {hamming_distance(p_genome, q_genome)}")
