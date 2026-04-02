#########################################################################################
# Genome Sequencing  - Bioinformatic II Course from Coursera
#
# Week 2 - String Reconstruction
# Input: An integer k followed by a list of k-mers Patterns.
# Output: A string Text with k-mer composition equal to Patterns.
# (If multiple answers exist, you may return any one.)
# Pseudocode:
# StringReconstruction(Patterns)
#     dB ← DeBruijn(Patterns)
#       Input format: ['GAGG', 'CAGG', 'GGGG', 'GGGA', 'CAGG', 'AGGG', 'GGAG']
#       Output format: {'GAG': ['AGG'], 'CAG': ['AGG', 'AGG'], 'GGG': ['GGG', 'GGA'], 'AGG': ['GGG'], 'GGA': ['GAG']}
#     path ← EulerianPath(dB)
#       Input format: ['0: 2', '1: 3', '2: 1', '3: 0 4', '6: 3 7', '7: 8', '8: 9', '9: 6']
#       Output format: 6->7->8->9->6->3->0->2->1->3->4
#     Text ← PathToGenome(path)
#       Input format: ['ACCGA', 'CCGAA', 'CGAAG', 'GAAGC', 'AAGCT']
#       Output format: ACCGAAGCT
#     return Text
# NOTE: YOU HAVE TO RUN THIS FROM THE BI2 FOLDER WITH THE COMMAND:
# python -m Wk2.Wk2_3_StringReconstruction (NO .py)
#########################################################################################

# Imports:
from Wk2.Wk2_2_EulerianPath import eulerian_path
from Wk1.Wk1_6_deBruijnGraph_kmer import debruijn_graph_from_kmers
from Wk1.Wk1_2_PathToGenome import path_to_genome

# String Reconstruction from patterns function
def string_reconstruction(kmers):
    """
    Recreates full genome from a string of kmers,
    where all the inputs are sanitized to be in the right format
    :param kmers -> kmers as a string separated by a space
    :return -> single genome
    """
    # Clean kmers for right format for dB
    kmers_sanitized = kmers.strip().split()
    dB = debruijn_graph_from_kmers(kmers_sanitized)

    # Clean output from dB
    dB_sanitized = [f"{key}: {' '.join(value)}" for key, value in dB.items()]
    path = eulerian_path(dB_sanitized)

    # Clean output from path
    path_sanitized = path.split("->")
    answer = path_to_genome([path_sanitized])

    return answer

if __name__ == "__main__":

# Sample dataset
    # kmers_orig = 'CTTA ACCA TACC GGCT GCTT TTAC'
    kmers = 'AAAT AATG ACCC ACGC ATAC ATCA ATGC CAAA CACC CATA CATC CCAG CCCA CGCT CTCA GCAT GCTC TACG TCAC TCAT TGCA'
    # Expected answer = CAAATGCATACGCTCATCACCCAG
    chullu = string_reconstruction(kmers)
    print(chullu)

# # From file
#     # Get dataset
#     from pathlib import Path as partho
#     current_dir = partho(__file__).parent
#     filename = input("Please enter the filename: ")
#     file_path = current_dir / filename

#     with open(file_path, 'r') as file:
#        # Read kmer length
#        k = int(file.readline().strip())
#        # Read the kmers
#        kmers = file.readline()

#     chullu = string_reconstruction(kmers)
#     print(chullu)


#     with open("Wk2/Wk2_3_output.txt", "w") as output_file:
#         output_file.write(chullu)