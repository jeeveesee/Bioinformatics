#########################################################################################
# Genome Sequencing  - Bioinformatic II Course from Coursera
#
# Week 2 - Generate contigs from a collection of reads (with imperfect coverage)
# Input: A collection of k-mers Patterns.
# Output: All contigs in DeBruijn(Patterns)
# Step 1: Create the debruijn graph from the kmers
# Step 2: Use maximal non-branching to find the kmers that form the contig
# Step 3: Create the contigs from the kmers
# NOTE: YOU HAVE TO RUN THIS FROM THE BI2 FOLDER WITH THE COMMAND:
# python -m Wk2.Wk2_9_Contigs_from_debruijn (NO .py)
#########################################################################################

# Imports:
from Wk2.Wk2_8_Maximal_Non_Branching_Paths import maximal_non_branching_paths
from Wk1.Wk1_6_deBruijnGraph_kmer import debruijn_graph_from_kmers
from Wk1.Wk1_2_PathToGenome import path_to_genome

# Generate contigs
def contigs(kmers):
    """
    Generates contigs from a list of kmer patterns
    It first finds the debruijn graph and then finds the maximal non-branching paths
    :param kmers -> kmers as a string separated by a space
    :return -> All contigs
    """
    # Step 1: Clean kmers for right format for dB and create dB graph
    kmers_sanitized = kmers.strip().split()
    dB = debruijn_graph_from_kmers(kmers_sanitized)

    # Step 2: Run maximnal non-branching path
    max_paths = maximal_non_branching_paths(dB)
    print(max_paths)

    #Step 3: Create the contigs
    answer = path_to_genome(max_paths)

    return answer

if __name__ == "__main__":

# # Sample dataset
#     kmers = 'ATG ATG TGT TGG CAT GGA GAT AGA'
#     bhaijo = contigs(kmers)
#     chullu = ' '.join(bhaijo)
#     print(chullu)

# From file
    # Get dataset
    from pathlib import Path as partho
    current_dir = partho(__file__).parent
    filename = input("Please enter the filename: ")
    file_path = current_dir / filename

    with open(file_path, 'r') as file:
       kmers = file.readline()

    bhaijo = contigs(kmers)
    chullu = ' '.join(bhaijo)
    # print(chullu)


    with open("Wk2/Wk2_9_output.txt", "w") as output_file:
        output_file.write(chullu)