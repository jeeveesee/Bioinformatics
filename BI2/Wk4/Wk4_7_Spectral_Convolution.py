#########################################################################################
# Genome Sequencing  - Bioinformatic II Course from Coursera
# Week 4 - Spectral Convolution problem
# Input: A collection of integers Spectrum in increasing order..
# Output: The list of elements in the convolution of Spectrum. If an element has multiplicity k,
# it should appear exactly k times; you may return the elements in any order.
# Hint: The convolution that you return should not contain any zeros
# A more efficient way might be to do the following:
"""
    n = len(spectrum)
    # Loop over each pair (i, j) where i > j to compute the differences
    for i in range(1, n):
        for j in range(i):
            diff = spectrum[i] - spectrum[j]
            if diff > 0:  # Only include positive differences (no zeros)
                convolution.append(diff)
"""
#########################################################################################

# Constants
# AMINO_ACID_MASSES = [i for i in range(57, 201)]

def spectral_convolution(experimental_spectrum):
    """
    Finds the spectral convolution of a given experimental spectrum

    Parameters:
        experimental_spectrum -> list of integers representing experimental spectrum masses

    Returns:
        list of amino acid masses in the convolutional spectrum
    """
    experimental_spectrum = sorted(experimental_spectrum)
    experimental_spectrum_conv = experimental_spectrum.copy()
    spectrum = []

    for mass in experimental_spectrum:
        for mass_conv in experimental_spectrum_conv:
            if (mass - mass_conv) > 0:
                spectrum.append(mass - mass_conv)
    return sorted(spectrum)


###########################################################################

if __name__ == "__main__":

# # Sample test
#     experimental_spectrum_raw = "0 137 186 323"
#     experimental_spectrum = list(map(int, experimental_spectrum_raw.split()))

#     # Expected answer = 137 137 186 186 323 49
#     conv = spectral_convolution(experimental_spectrum)
#     # formatted = ["-".join(map(str, p)) for p in leaders]
#     formatted = " ".join(map(str, conv))
#     # print("Leader peptides are:\n")
#     # print("\n".join(formatted))
#     print(formatted)
#     # print("\nBest score is:", best_score)
#     # print(*answer)

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
#         # N = int(file.readline())
#         experimental_spectrum_raw = file.readline()
#     experimental_spectrum = list(map(int, experimental_spectrum_raw.split()))

#     conv = spectral_convolution(experimental_spectrum)
#     answer = " ".join(map(str, conv))
#     # print(answer)

#     with open("Wk4_7_output.txt", "w") as output_file:
#         # output_file.write(str(answer))
#         output_file.write(answer)


# # For Exercise Break

#     # # Get dataset
#     # from pathlib import Path as partho
#     # current_dir = partho(__file__).parent
#     # filename = input("Please enter the filename: ")
#     # file_path = current_dir / filename

#     # #NOTE: For Wk3_2 exercise break, we need to concatenate all the kmers on different lines
#     # with open(file_path, 'r') as file:
#     #    dna_string = file.read().replace('\n', '')
#     # # print(dna_string)

#     # For Wk4_4  Cyclopeptide Leaderboard
#     # Spectrum(25) for Tyrocidine B1
#     experimental_spectrum_raw = "0 97 99 113 114 115 128 128 147 147 163 186 227 241 242 244 244 256 260 261 262 283 291 309 330 333 340 347 385 388 389 390 390 405 435 447 485 487 503 504 518 544 552 575 577 584 599 608 631 632 650 651 653 672 690 691 717 738 745 770 779 804 818 819 827 835 837 875 892 892 917 932 932 933 934 965 982 989 1039 1060 1062 1078 1080 1081 1095 1136 1159 1175 1175 1194 1194 1208 1209 1223 1322"
#     experimental_spectrum = list(map(int, experimental_spectrum_raw.split()))
#     N = 1000

#     answer = leaderboard_cyclopeptide_sequencing(experimental_spectrum, N)
#     print('-'.join(answer))

#     # with open("Wk3_2_exercisebreak_output.txt", "w") as output_file:
#     #     output_file.write('\n'.join(answer))

# EXAM
    experimental_spectrum_raw = "0 86 160 234 308 320 382"
    experimental_spectrum = list(map(int, experimental_spectrum_raw.split()))

    # Expected answer = 137 137 186 186 323 49
    conv = spectral_convolution(experimental_spectrum)
    formatted = " ".join(map(str, conv))
    print(formatted)