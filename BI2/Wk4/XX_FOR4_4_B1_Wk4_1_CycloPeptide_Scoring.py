#########################################################################################
# Genome Sequencing  - Bioinformatic II Course from Coursera
#
# Week 4 - Cyclopeptide scoring - noisy experimental spectrum
# Input: An amino acid string Peptide and a collection of integers Spectrum
# Output: The score of Peptide against Spectrum, Score(Peptide, Spectrum)
# The scoring function should take into account the multiplicities of shared masses
# In general, if a mass occurs m times in the theoretical spectrum of Peptide
# and n times in the experimental spectrum Spectrum,
# then it contributes the minimum of m and n to Score(Peptide, Spectrum).
# NOTE: YOU HAVE TO RUN THIS FROM THE BI2 FOLDER WITH THE COMMAND:
# python -m Wk4.Wk4_1_CycloPeptide_Scoring (NO .py)
#########################################################################################

# Use mass-based cyclospectrum from Wk3.Wk3_8_CycloPeptide_Sequencing
from Wk3.Wk3_8_CycloPeptide_Sequencing import cyclospectrum


# Cyclopeptide Scoring (expects peptide as a tuple of integer masses)
def cyclopeptide_scoring(peptide: tuple[int, ...], experimental_spectrum: list[int]) -> int:
    """
    Compute the cyclic score of a peptide represented as a tuple of masses.

    Parameters:
    peptide -> tuple of integer masses
    experimental_spectrum -> list of integers representing experimental spectrum

    Returns:
    integer score (counts multiplicities)
    """
    from collections import Counter

    theoretical_spectrum = cyclospectrum(peptide)

    theo_count = Counter(theoretical_spectrum)
    exper_count = Counter(experimental_spectrum)

    score = 0
    for mass in theo_count:
        if mass in exper_count:
            score += min(theo_count[mass], exper_count[mass])

    return score


###########################################################################

if __name__ == "__main__":

# Sample test
    peptide_string = 'NQEL'
    experimental_spectrum = [0, 99, 113, 114, 128, 227, 257, 299, 355, 356, 370, 371, 484]
    # Expected answer = 11
    answer = cyclopeptide_scoring(peptide_string, experimental_spectrum)
    # print('\n'.join(answer))
    print(answer)

# # From file

#     # Get dataset
#     from pathlib import Path as partho
#     current_dir = partho(__file__).parent
#     filename = input("Please enter the filename: ")
#     file_path = current_dir / filename

#     #NOTE: For Wk3_2, there is a newline character after the peptide string, so remove that
#     with open(file_path, 'r') as file:
#         # Remove trailing newline/whitespace from the peptide string (was causing '\n' to be iterated)
#         peptide_string = file.readline().strip()
#         experimental_spectrum_raw = file.readline()

#     experimental_spectrum = list(map(int, experimental_spectrum_raw.split()))

#     answer = cyclopeptide_scoring(peptide_string, experimental_spectrum)
#     # print('\n'.join(answer))
#     print(answer)

#     with open("Wk4_1_output.txt", "w") as output_file:
#         output_file.write(str(answer))


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

#     answer = peptide_encoding(dna_string, peptide_string)
#     # print(len(answer))

#     with open("Wk3_2_exercisebreak_output.txt", "w") as output_file:
#         output_file.write('\n'.join(answer))
#         # Answer is 0!!