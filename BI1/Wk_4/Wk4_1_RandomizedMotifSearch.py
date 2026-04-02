#########################################################################################
# TASK:
# Use random motif search method to generate  profiles to
# find a set of motifs that are likely implanted by a common motif
# NOTE: YOU HAVE TO RUN THIS FROM THE BI1 FOLDER WITH THE COMMAND:
# python -m Wk_4.Wk_4_1_RandomizedMotifSearch
##########################################################################################

# def randomized_motif_search(dna, k, t):

import random
from Wk_3.Wk3_7_greedy_motif_search_w_pseudocounts_EFFICIENT import most_probable_kmer, create_profile_matrix, score_motifs

def randomized_motif_search(dna_list, k, t):
    # Step 1: Randomly select k-mers from each string in Dna to form initial motifs
    random_motifs = []
    for dna in dna_list:
        start_index = random.randint(0, len(dna) -k)
        random_motif = dna[start_index:start_index + k]
        random_motifs.append(random_motif)

    best_motifs = random_motifs
    best_score = score_motifs(best_motifs)

    while True:
        current_profile = create_profile_matrix(best_motifs)
        new_motifs = []

        for dna in dna_list:
            _, next_motif = most_probable_kmer(dna, k, current_profile)
            new_motifs.append(next_motif)

        current_score = score_motifs(new_motifs)

        if current_score < best_score:
            best_score = current_score
            best_motifs = new_motifs
        else:
            return best_motifs


if __name__ == "__main__":
    # Example from text
    # dna_list = ['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA', 'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG', 'TAGTACCGAGACCGAAAGAAGTATACAGGCGT', 'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC', 'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']
    # k = 8
    # t = 5    # print(randomized_motif_search(dna_list, k, t))

   # Get dataset
   from pathlib import Path as partho
   current_dir = partho(__file__).parent
   filename = input("Please enter the filename: ")
   file_path = current_dir / filename

   with open(file_path, 'r') as file:
       # Read k and t (Row 0)
       k_and_t = file.readline().strip()
       k, t = map(int, k_and_t.split())

       # Read the DNA strings
       dna_list_raw = file.readline().strip()
       dna_list = dna_list_raw.split()
       #print(dna_list)
   all_best_motifs = []

   for i in range(1000):
       best_motifs = randomized_motif_search(dna_list, k, t)
       all_best_motifs.append(best_motifs)
       if (i+1) % 100 == 0:
          print(f"Iteration {i+1}")
    
   best_answer = min(all_best_motifs, key=score_motifs)
   print("\n\nThe best motifs are: ", *best_answer)