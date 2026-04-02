#########################################################################################
# Genome Sequencing  - Bioinformatic II Course from Coursera
#
# Week 1 - Overlap Graph problem
# Input: A collection Patternsof k-mers
# Output: The overlap graph Overlap(Patterns), in the form of an adjacency list. 
# (You may return the nodes and their edges in any order.)
# NOTE: Now the input kmers are NOT in the correct order
#########################################################################################

# Overlap Graph function
##################################

def overlap_graph(kmer_list):
    overlap_dict = {}

    for kmer1 in kmer_list:
        for kmer2 in kmer_list:
            if kmer1 != kmer2: # Avoiding self loops
                # Overlap of k-1 suffix of kmer1 = k-1 prefix of kmer2
                if kmer1[1:] == kmer2[:-1]:
                    if kmer1 in overlap_dict:
                        overlap_dict[kmer1].append(kmer2)
                    else:
                        overlap_dict[kmer1] = [kmer2] # Add kmer2 as a list, so we can append to it later
    return overlap_dict


# Get dataset
##################################
if __name__ == '__main__':

    # Sample input
    kmer_list = ['ATGCG', 'GCATG', 'CATGC', 'AGGCA', 'GGCAT', 'GGCAC']
    print(overlap_graph(kmer_list))
    answer = overlap_graph(kmer_list)
    for key, values in answer.items():
        # Join the list of values into a single string with ', ' as the separator
        values_str = ", ".join(values)
        # Print the key, a colon and space, and the formatted values
        print(f"{key}: {values_str}")


    # # Get dataset
    # from pathlib import Path as partho
    # current_dir = partho(__file__).parent
    # filename = input("Please enter the filename: ")
    # file_path = current_dir / filename

    # with open(file_path, 'r') as file:
    #    # Read the kmer strings
    #    kmer_list_raw = file.readline().strip()
    #    kmer_list = kmer_list_raw.split()

    # answer = overlap_graph(kmer_list)

    # with open("Wk1_3_output.txt", "w") as output_file:
    #     for key, values in answer.items():

    #         values_str = ", ".join(values)
    #         line = f"{key}: {values_str}\n"
    #         output_file.write(line)






