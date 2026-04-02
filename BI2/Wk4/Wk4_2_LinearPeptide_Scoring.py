#########################################################################################
# Genome Sequencing  - Bioinformatic II Course from Coursera
#
# Week 4 - Linear peptide Scoring - noisy experimental spectrum
# Input: An amino acid string Peptide and a collection of integers Spectrum.
# Output: The linear score of Peptide with respect to Spectrum, LinearScore(Peptide, Spectrum).
# Step 1: Create LinearScore(Peptide, Spectrum) that computes the linear score of Peptide
# This is part 1 of the Trim Leaderboard algorithm
# NOTE: YOU HAVE TO RUN THIS FROM THE BI2 FOLDER WITH THE COMMAND:
# python -m Wk4.Wk4_2_LinearPeptide_Scoring (NO .py)
#########################################################################################

# Imports:
from Wk3.Wk3_4_LinearSpectrum import linear_spectrum

# Constants
AMINO_ACID_MASS_MAP = {
    'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101,
    'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115, 'K': 128,
    'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156,
    'Y': 163, 'W': 186}

# Linear Scoring
def linear_scoring(peptide_string: str, experimental_spectrum: list[int]):
    """
    Takes an peptide sequence and experimental spectrum to find its linear score

    Parameters:
    peptide_string -> string of amino acid letters
    experimental_spectrum -> list of integers representing experimental spectrum

    Returns:
    integer score of linear against experimental spectrum (takes multiplicities into account)
    """

    from collections import Counter

    # Get theoretical linear spectrum
    theoretical_spectrum = linear_spectrum(peptide_string)

    # Count occurrences in both spectra
    theo_count = Counter(theoretical_spectrum)
    exper_count = Counter(experimental_spectrum)

    # Calculate score based on minimum occurrences
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
    # Expected answer = 8
    # EMPTY: 0 - YES
    # N: 114 - YES
    # Q: 128 - YES
    # E: 129 - NOT PRESENT
    # L: 113 - YES
    # NQ: 114 + 128 = 242 - NOT PRESENT
    # QE: 128 + 129 = 258 - YES
    # EL: 129 + 113 = 242 - NOT PRESENT
    # NQE: 114 + 128 + 129 = 371 - YES
    # QEL: 128 + 129 + 113 = 370 - YES
    # NQEL: 114 + 128 + 129 + 113 = 484 - YES
    answer = linear_scoring(peptide_string, experimental_spectrum)
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

#     answer = linear_scoring(peptide_string, experimental_spectrum)
#     # print('\n'.join(answer))
#     print(answer)

#     with open("Wk4_2_output.txt", "w") as output_file:
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

#     answer = peptide_encoding(RNA_CODON_MAP, dna_string, peptide_string)
#     # print(len(answer))

#     with open("Wk3_2_exercisebreak_output.txt", "w") as output_file:
#         output_file.write('\n'.join(answer))
#         # Answer is 0!!

# # EXAM
#     peptide_string = 'PEEP'
#     experimental_spectrum = [0, 97, 129, 129, 129, 194, 226, 323, 323, 355, 452]
#     answer = linear_scoring(peptide_string, experimental_spectrum)
#     # print('\n'.join(answer))
#     print(answer)