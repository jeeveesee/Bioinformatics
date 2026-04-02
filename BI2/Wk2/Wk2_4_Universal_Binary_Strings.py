#########################################################################################
# Genome Sequencing  - Bioinformatic II Course from Coursera
#
# Week 2 - Universal Circular Binary Strings
#   Input: An integer k
#   Output: A k-universal circular string
# Pseudocode:
# Let BinaryStringsk be the set of all 2^k binary k-mers.
# The only thing we need to do is solve the k-Universal Circular String Problem
# is to find an Eulerian cycle in DeBruijn(BinaryStringsk).
# NOTE: YOU HAVE TO RUN THIS FROM THE BI2 FOLDER WITH THE COMMAND:
# python -m Wk2.Wk_2_4_Universal_Binary_Strings (NO .py)
#########################################################################################
#########################################################################################

# Imports
from itertools import product
from Wk1.Wk1_6_deBruijnGraph_kmer import debruijn_graph_from_kmers
from Wk2.Wk2_2_EulerianPath import eulerian_path
from Wk1.Wk1_2_PathToGenome import path_to_genome

# Create set of 2^k binary k-mers
def binary_kmers(k):
    for tup in product('01', repeat=k):
        yield ''.join(tup)

def universal_binary_string(k):
    """
    Recreates full genome from a string of kmers,
    where all the inputs are sanitized to be in the right format
    :param kmers -> kmers as a string separated by a space
    :return -> single genome
    """
    # Create binary kmers
    kmers = list(set(binary_kmers(k)))

    # Clean kmers for right format for dB
    # kmers_sanitized = kmers.split().strip()
    dB = debruijn_graph_from_kmers(kmers)

    # Clean output from dB
    dB_sanitized = [f"{key}: {' '.join(value)}" for key, value in dB.items()]
    path = eulerian_path(dB_sanitized)

    # Clean output from path
    path_sanitized = path.split("->")
    answer = path_to_genome(path_sanitized)

    return answer


if __name__ == "__main__":

# # Sample dataset
    # k = 3
    # answer = universal_binary_string(k)
    # print(answer[:-k+1])

# # From file
    # Get dataset
    from pathlib import Path as partho
    current_dir = partho(__file__).parent
    filename = input("Please enter the filename: ")
    file_path = current_dir / filename

    with open(file_path, 'r') as file:
       # Read kmer length
       k = int(file.readline().strip())

    chullu = universal_binary_string(k)
    print(chullu[:-k+1])
    # The above is done to close the loop and make it circular
    # as the overlap is k-1, so you leave out k-1 nucleotides at the end


    with open("Wk2/Wk2_4_output.txt", "w") as output_file:
        output_file.write(chullu[:-k+1])




