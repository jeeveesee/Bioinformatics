#########################################################################################
# Genome Sequencing  - Bioinformatic II Course from Coursera
#
# Week 1 - de Bruijn Graphs
# Input: A collection of k-mers Patterns
# Output: The adjacency list of the de Bruijn graph DeBruijn(Patterns)
# NOTE:
# This is NOT the same as the Overlap Graph from a list of kmers!
# It is also NOT the same as the de Bruijn graph from a single string of DNA
# This is from a list of kmers!!! DIFFERENT FROM Wk1_5_deBruijnGraph.py
#########################################################################################

# de Bruijn Graph function
##################################

def debruijn_graph_from_kmers(kmers):
    """
    Creates a de Bruijn Graph from list of kmers of same size
    :param kmers -> list of string of same size kmers
    :return -> adjacency list from graph in a dictionary
    """
    debruijn_edges = {}

    for kmer in kmers:
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
    kmers = ['GAGG', 'CAGG', 'GGGG', 'GGGA', 'CAGG', 'AGGG', 'GGAG']
    answer = debruijn_graph_from_kmers(kmers)
    print(answer)
    for key, values in sorted(answer.items()):
        values_str = " ".join(sorted(values))
        print(f"{key}: {values_str}")


    # # Get dataset
    # from pathlib import Path as partho
    # current_dir = partho(__file__).parent
    # filename = input("Please enter the filename: ")
    # file_path = current_dir / filename

    # with open(file_path, 'r') as file:
    # #    # Read the kmer integer
    # #    k = int(file.readline().strip())
    #    # Read the DNA string
    #    kmers = file.readline().strip().split()

    # answer = debruijn_graph_from_kmers(kmers)

    # with open("Wk1_6_output.txt", "w") as output_file:
    #     for key, values in sorted(answer.items()):

    #         values_str = " ".join(sorted(values))
    #         line = f"{key}: {values_str}\n"
    #         output_file.write(line)