#########################################################################################
# Genome Sequencing  - Bioinformatic II Course from Coursera
#
# Week 1 - String Composition problem
# For a given DNA read and an integer of length k,
# generate all the kmers of length k
#########################################################################################

def string_composition(dna_read, kmer_len):
    kmers = []
    for i in range(len(dna_read) - kmer_len + 1):
        kmers.append(dna_read[i:i+kmer_len])
    return kmers

if __name__ == '__main__':
    # Sample input
    # dna_read = "CAATCCAAC"
    # kmer_len = 5
    # print(string_composition(dna_read, kmer_len))

    # Get dataset
    from pathlib import Path as partho
    current_dir = partho(__file__).parent
    filename = input("Please enter the filename: ")
    file_path = current_dir / filename

    with open(file_path, 'r') as file:
       # Read k (Row 0)
       k = file.readline().strip()
       # Read the DNA strings
       dna_list_raw = file.readline().strip()
       dna_list = dna_list_raw.split()

    answer = string_composition(dna_list[0], int(k))
    answer_string = ' '.join(answer)

    with open("Wk1_1_output.txt", "w") as output_file:
        output_file.writelines(answer_string)