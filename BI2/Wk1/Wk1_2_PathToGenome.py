#########################################################################################
# Genome Sequencing  - Bioinformatic II Course from Coursera
#
# Week 1 - Path to Genome problem
# Input: A sequence path of k-mers Pattern1, … ,Patternn
# such that the last k - 1 symbols of Patterni are equal
# to the first k-1 symbolsof Patterni+1 for 1 ≤ i ≤ n-1
# Output: A string Text of length k+n-1 such that the
# i-th # k-mer in Text is equal to Patterni (for 1 ≤ i ≤ n).
# NOTE: We're assuming that the input kmers are in the correct order
#########################################################################################

# Path to Genome function
##################################
def path_to_genome(list_of_kmer_list):
    """
    Creates path to a genome through a list of kmers
    :param kmer_list -> list of string of same size kmers in a path
    :return -> reconstructed genome with overlap of 1
    """
    genomes = []
    for kmer_list in list_of_kmer_list:
        genome = kmer_list[0] # Starting with the first kmer in full
        kmer_len = len(kmer_list[0])
        num_patterns = len(kmer_list)

        for i in range(1, num_patterns):
            genome += kmer_list[i][(kmer_len - 1):]
        genomes.append(genome)

    return genomes

# Get dataset
##################################
if __name__ == '__main__':
    # Sample input
    list_of_kmer_list = [['ACCGA', 'CCGAA', 'CGAAG', 'GAAGC', 'AAGCT'], ['ACCGA', 'CCGAA', 'CGAAG', 'GAAGC', 'AAGCT']]
    print(path_to_genome(list_of_kmer_list))


    # # Get dataset
    # from pathlib import Path as partho
    # current_dir = partho(__file__).parent
    # filename = input("Please enter the filename: ")
    # file_path = current_dir / filename

    # with open(file_path, 'r') as file:
    #    # Read the kmer strings
    #    kmer_list_raw = file.readline().strip()
    #    kmer_list = kmer_list_raw.split()

    # answer = path_to_genome(kmer_list)
    # #answer_string = ' '.join(answer)

    # with open("Wk1_2_output.txt", "w") as output_file:
    #     output_file.writelines(answer)

