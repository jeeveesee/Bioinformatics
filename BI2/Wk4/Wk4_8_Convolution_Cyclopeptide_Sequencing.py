#########################################################################################
# Genome Sequencing  - Bioinformatic II Course from Coursera
#
# Week 4 - Convolution Cyclopeptide Sequencing
# Input: An integer M, an integer N, and a collection of (possibly repeated) integers Spectrum.
# Output: A cyclic peptide LeaderPeptide with amino acids taken only from the
# top M elements (and ties) of the convolution of Spectrum that fall between
# 57 and 200, and where the size of Leaderboard is restricted to the top N (and ties).

"""
We now have the outline for a new cyclopeptide sequencing algorithm.
Given an experimental spectrum, we first compute the convolution of an experimental spectrum.
We then select the M most frequent elements between 57 and 200 in the convolution
to form an extended alphabet of candidate amino acid masses.
In order to be fair, we should include the top M elements of the convolution "with ties".
Finally, we run the algorithm LeaderboardCyclopeptideSequencing,
where the amino acid masses are restricted to this alphabet.
We call this algorithm ConvolutionCyclopeptideSequencing.
"""

# NOTE: YOU HAVE TO RUN THIS FROM THE BI2 FOLDER WITH THE COMMAND:
# python -m Wk4.Wk4_8_Convolution_Cyclopeptide_Sequencing (NO .py)
#########################################################################################

# Imports:
from Wk4.Wk4_6_Leaderboard_Cyclopeptide_Sequencing_Extended import leaderboard_cyclopeptide_sequencing
from Wk4.Wk4_7_Spectral_Convolution import spectral_convolution

def convolution_cyclopeptide_sequencing(experimental_spectrum, M, N):
    """
    Finds the convolution cyclopeptide sequencing of a given experimental spectrum

    Parameters:
        experimental_spectrum -> list of integers representing experimental spectrum masses
        M -> integer, number of top elements to consider from convolution
        N -> integer, size of leaderboard

    Returns:
        list of leader peptides with amino acid masses in the convolutional spectrum
    """
    # Get convolution of the experimental spectrum
    conv_spectrum = spectral_convolution(experimental_spectrum)

    # Count frequencies of each mass in the convolution
    from collections import Counter
    mass_counts = Counter([mass for mass in conv_spectrum if 57 <= mass <= 200])
    # print(f"{mass_counts=}")

    # Get the M most common masses (with ties)
    if not mass_counts:
        return [], 0

    most_common = mass_counts.most_common()
    # print(f"\n{most_common=}")
    # What is the smallest maximum count of the masses that we should care about when placed in descending order
    threshold_count = most_common[M-1][1] if M <= len(most_common) else 1

    extended_alphabet = [mass for mass, count in most_common if count >= threshold_count]

    # Run leaderboard cyclopeptide sequencing with the extended alphabet
    leaders, best_score = leaderboard_cyclopeptide_sequencing(experimental_spectrum, N, extended_alphabet)

    return leaders, best_score




###########################################################################

if __name__ == "__main__":

# Sample test
    experimental_spectrum_raw = "57 57 71 99 129 137 170 186 194 208 228 265 285 299 307 323 356 364 394 422 493"
    experimental_spectrum = list(map(int, experimental_spectrum_raw.split()))
    M = 20
    N = 60
    # Expected answer = 99-71-137-57-72-57
    leaders, best_score = convolution_cyclopeptide_sequencing(experimental_spectrum, M, N)
    formatted = ["-".join(map(str, p)) for p in leaders]
    print("Leader peptides are:\n")
    print("\n".join(formatted))
    print("\nBest score is:", best_score)
    # print(*answer)

# # From file

#     # Get dataset
#     from pathlib import Path as partho
#     current_dir = partho(__file__).parent
#     filename = input("Please enter the filename: ")
#     file_path = current_dir / filename

#     #NOTE: For Wk3_2, there is a newline character after the peptide string, so remove that
#     with open(file_path, 'r') as file:
#         # Read peptide string and remove ALL whitespace (spaces/newlines) — dataset may contain spaces
#         # leaderboard = file.readline().split()
#         M = int(file.readline())
#         N = int(file.readline())
#         experimental_spectrum_raw = file.readline()
#     experimental_spectrum = list(map(int, experimental_spectrum_raw.split()))

#     leaders, best_score = convolution_cyclopeptide_sequencing(experimental_spectrum, M, N)
#     formatted = ["-".join(map(str, p)) for p in leaders]
#     answer = "\n".join(formatted)

#     with open("Wk4/Wk4_8_output.txt", "w") as output_file:
#         # output_file.write(str(answer))
#         output_file.write(answer)


# # # For Exercise Break
#     # For Wk4_8  Convolution  Cyclopeptide Sequencing
#     M = 20
#     N = 1000
#     experimental_spectrum_raw = "0 97 99 113 114 115 128 128 147 147 163 186 227 241 242 244 244 256 260 261 262 283 291 309 330 333 340 347 385 388 389 390 390 405 435 447 485 487 503 504 518 544 552 575 577 584 599 608 631 632 650 651 653 672 690 691 717 738 745 770 779 804 818 819 827 835 837 875 892 892 917 932 932 933 934 965 982 989 1039 1060 1062 1078 1080 1081 1095 1136 1159 1175 1175 1194 1194 1208 1209 1223 1322"
#     experimental_spectrum = list(map(int, experimental_spectrum_raw.split()))

#     leaders, best_score = convolution_cyclopeptide_sequencing(experimental_spectrum, M, N)
#     formatted = ["-".join(map(str, p)) for p in leaders]
#     answer = " ".join(formatted)
#     print(answer)
