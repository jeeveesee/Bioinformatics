#########################################################################################
# Genome Sequencing  - Bioinformatic III Course from Coursera
#
# Week 5 - Shared k-mers
# Shared k-mers Problem: Given two strings, find all their shared k-mers.
# Use Wk5_9 for BIG BACTERIA GENOMES!!
#########################################################################################

from pathlib import Path as partho
from collections import defaultdict


def reverse_complement(dna):
    """Return the reverse complement of a DNA string.

    Parameters:
        dna (str): A DNA string composed of A, T, C, G characters.

    Returns:
        str: The reverse complement of the input DNA string.
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(dna))


def build_kmer_index(dna, k):
    """Build a dictionary mapping each k-mer in a DNA string to its starting positions.

    Parameters:
        dna (str): A DNA string.
        k (int): The length of k-mers to index.

    Returns:
        dict: A mapping from k-mer string to a list of integer starting positions.
    """
    index = defaultdict(list)
    for i in range(len(dna) - k + 1):
        kmer = dna[i:i + k]
        index[kmer].append(i)
    return index


def shared_kmers(k, dna1, dna2):
    """Find all shared k-mers between two DNA strings, including reverse complement matches.

    Two k-mers are considered shared if one equals the other or its reverse complement.
    For each starting position i in dna1, this function first records positions where
    the exact k-mer appears in dna2, then positions where the reverse complement appears.

    Parameters:
        k (int): The length of k-mers to compare.
        dna1 (str): The first DNA string.
        dna2 (str): The second DNA string.

    Returns:
        str: A newline-separated string of (i, j) pairs where dna1[i:i+k] equals
             dna2[j:j+k] or its reverse complement.
    """
    index2 = build_kmer_index(dna2, k)
    pairs = []

    for i in range(len(dna1) - k + 1):
        kmer = dna1[i:i + k]
        # Direct matches in dna2
        for j in index2.get(kmer, []):
            pairs.append((i, j))
        # Reverse complement matches in dna2
        rc = reverse_complement(kmer)
        if rc != kmer:  # avoid duplicating palindromic k-mers
            for j in index2.get(rc, []):
                pairs.append((i, j))

    return '\n'.join(f"({i}, {j})" for i, j in pairs)


###########################################################################

if __name__ == "__main__":
    # # Sample test
    # kmer = 3
    # dna_string_1 = "AAACTCATC"
    # dna_string_2 = "TTTCAAATC"

    # # Expected answer =
    # # (0, 4)
    # # (0, 0)
    # # (4, 2)
    # # (6, 6)
    # answer = shared_kmers(kmer, dna_string_1, dna_string_2)
    # print(answer)

    # From file
    current_dir_1 = partho(__file__).parent
    filename_1 = input("Please enter the filename: ")
    file_path_1 = current_dir_1 / filename_1

    with open(file_path_1, "r") as file:
        dna_string_1 = file.readline().strip()

    current_dir_2 = partho(__file__).parent
    filename_2 = input("Please enter the filename: ")
    file_path_2 = current_dir_2 / filename_2

    with open(file_path_2, "r") as file:
        dna_string_2 = file.readline().strip()

    k_val = 30

    answer = shared_kmers(k_val, dna_string_1, dna_string_2)
    # print(answer)

    current_dir = partho(__file__).parent
    with open(current_dir / "Wk5_10_BIGoutput.txt", "w") as output_file:
        output_file.write(answer)

"""
 ---
  Walkthrough: k=3, dna1="AAACTCATC", dna2="TTTCAAATC"

  Step 1 — Index all 3-mers in dna2:

  ┌──────────┬───────┐
  │ Position │ k-mer │
  ├──────────┼───────┤
  │ 0        │ TTT   │
  ├──────────┼───────┤
  │ 1        │ TTC   │
  ├──────────┼───────┤
  │ 2        │ TCA   │
  ├──────────┼───────┤
  │ 3        │ CAA   │
  ├──────────┼───────┤
  │ 4        │ AAA   │
  ├──────────┼───────┤
  │ 5        │ AAT   │
  ├──────────┼───────┤
  │ 6        │ ATC   │
  └──────────┴───────┘

  Step 2 — For each k-mer in dna1, check for direct matches AND reverse complement matches in the index:


  ┌─────┬─────────────┬───────────────────────┬─────────┬────────────────────────┬──────────────┐
  │  i  │ dna1[i:i+3] │ Direct match in dna2? │ RevComp │ RevComp match in dna2? │ Pairs added  │
  ├─────┼─────────────┼───────────────────────┼─────────┼────────────────────────┼──────────────┤
  │ 0   │ AAA         │ pos 4                 │ TTT     │ pos 0                  │ (0,4), (0,0) │
  ├─────┼─────────────┼───────────────────────┼─────────┼────────────────────────┼──────────────┤
  │ 1   │ AAC         │ —                     │ GTT     │ —                      │ —            │
  ├─────┼─────────────┼───────────────────────┼─────────┼────────────────────────┼──────────────┤
  │ 2   │ ACT         │ —                     │ AGT     │ —                      │ —            │
  ├─────┼─────────────┼───────────────────────┼─────────┼────────────────────────┼──────────────┤
  │ 3   │ CTC         │ —                     │ GAG     │ —                      │ —            │
  ├─────┼─────────────┼───────────────────────┼─────────┼────────────────────────┼──────────────┤
  │ 4   │ TCA         │ pos 2                 │ TGA     │ —                      │ (4,2)        │
  ├─────┼─────────────┼───────────────────────┼─────────┼────────────────────────┼──────────────┤
  │ 5   │ CAT         │ —                     │ ATG     │ —                      │ —            │
 ─┼───────────────────────┼─────────┼────────────────────────┼──────────────┤
  │ 6   │ ATC         │ pos 6                 │ GAT     │ —                      │ (6,6)        │

"""