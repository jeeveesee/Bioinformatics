#########################################################################################
# Genome Sequencing  - Bioinformatic II Course from Coursera
#
# Week 1 - de Bruijn Graphs
# Input: An integer k and a string Text.
# Output: DeBruijnk(Text), in the form of an adjacency list.
# NOTE:
# This is NOT the same as the Overlap Graph from a list of kmers!
# This is from a single string of DNA

# NOTE:
# Hamiltonian → nodes,
# Eulerian → edges

#########################################################################################

# de Bruijn Graph function
##################################

def debruijn_graph(k, dna):
    debruijn_edges = {}

    for i in range(len(dna) - k + 1):
        kmer = dna[i:i+k]
        prefix = kmer[:-1]
        suffix = kmer[1:]

        if prefix in debruijn_edges:
            debruijn_edges[prefix].append(suffix)
        else:
            debruijn_edges[prefix] = [suffix] # Add suffix as a list, so we can append to it later
    return debruijn_edges


# Get dataset
##################################
if __name__ == '__main__':

    # Sample input
    k = 3
    dna = 'TAATGGGATGCCATGTT'
    answer = debruijn_graph(k, dna)
    for key, values in answer.items():
        values_str = " ".join(values)
        print(f"{key}: {values_str}")


    # # Get dataset
    # from pathlib import Path as partho
    # current_dir = partho(__file__).parent
    # filename = input("Please enter the filename: ")
    # file_path = current_dir / filename

    # with open(file_path, 'r') as file:
    #    # Read the kmer integer
    #    k = int(file.readline().strip())
    #    # Read the DNA string
    #    dna = file.readline().strip()

    # answer = debruijn_graph(k, dna)

    # with open("Wk1_5_output.txt", "w") as output_file:
    #     for key, values in sorted(answer.items()):

    #         values_str = " ".join(values)
    #         line = f"{key}: {values_str}\n"
    #         output_file.write(line)