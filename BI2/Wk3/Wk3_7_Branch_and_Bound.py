#########################################################################################
# Genome Sequencing  - Bioinformatic II Course from Coursera
#
# Week 3 - branch and bound LINEAR peptides
# Input: An integer n.
# Output: The number of subpeptides of a linear peptide of length n

"""
More generally, brute force algorithms that enumerate all candidate solutions
but discard large subsets of hopeless candidates by using various consistency
conditions are known as branch-and-bound algorithms.
Each such algorithm consists of a branching step to increase the
number of candidate solutions, followed by a bounding step to remove hopeless candidates.
In our branch-and-bound algorithm for the Cyclopeptide Sequencing Problem,
the branching step will extend each candidate peptide of length k
into 18 peptides of length k + 1, and the bounding step will remove inconsistent peptides
from consideration.

Note that the spectrum of a linear peptide does not contain as many masses as the spectrum
of a cyclic peptide with the same amino acid sequence.
For instance, the theoretical spectrum of the cyclic peptide
NQEL contains 14 masses
(corresponding to "", N, Q, E, L, LN, NQ, QE, EL, ELN, LNQ, NQE, QEL, and NQEL).
However, the theoretical spectrum of the linear peptide NQEL, shown below,
does not contain masses corresponding to
LN, LNQ, or ELN, since these subpeptides “wrap around” the end of the linear peptide.
"""
# The answer can be calculated as [n*(n+1)/2]+1
"""
For a linear peptide with n amino acids:

Length 1 subpeptides: n choices (positions 1, 2, 3, ..., n)
Length 2 subpeptides: n-1 choices (positions 1-2, 2-3, ..., (n-1)-n)
Length 3 subpeptides: n-2 choices (positions 1-2-3, 2-3-4, ..., (n-2)-(n-1)-n)
...
Length n subpeptides: 1 choice (the entire peptide) (n-n or 1 full choice)
Empty peptide: 1 choice

Total = 1 + n + (n-1) + (n-2) + ... + 2 + 1 = 1 + (sum of 1st n numbers)
= 1 + (n*(n-1)/2)
"""
#########################################################################################

def count_subpeptides_LINEAR(n):
    return ((n*(n+1)/2))+1


###########################################################################

if __name__ == "__main__":

# # Sample test
#     n = 4
#     # Expected answer = 11
#     answer = count_subpeptides_LINEAR(n)
#     ## print(''.join(answer))
#     print(answer)


# # From file

    # Get dataset
    from pathlib import Path as partho
    current_dir = partho(__file__).parent
    filename = input("Please enter the filename: ")
    file_path = current_dir / filename

    #NOTE: Make sure to remove the extra line at the end of the file
    with open(file_path, 'r') as file:
       n = file.readline()

    answer = count_subpeptides_LINEAR(int(n))
    # print('\n'.join(answer))

    with open("Wk3_7_output.txt", "w") as output_file:
    #   output_file.write('\n'.join(answer))
    #   output_file.write(" ".join(map(str, answer)))
        output_file.write(str(answer))



# # For Exercise Break

#     # Get dataset
#     from pathlib import Path as partho
#     current_dir = partho(__file__).parent
#     filename = input("Please enter the filename: ")
#     file_path = current_dir / filename

#     #NOTE: For Wk3_2 exercise break, we need to concatenate all the kmers on different lines
#     with open(file_path, 'r') as file:
#        dna_string = file.read().replace('\n', '')
#     # print(dna_string)

#     # Val-Lys-Leu-Phe-Pro-Trp-Phe-Asn-Gln-Tyr
#     peptide_string = 'VKLFPWFNQY'

#     answer = peptide_encoding(RNA_CODON_MAP, dna_string, peptide_string)
#     # print(len(answer))

#     with open("Wk3_2_exercisebreak_output.txt", "w") as output_file:
#         output_file.write('\n'.join(answer))
#         # Answer is 0!!