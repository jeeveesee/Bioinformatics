#########################################################################################
# Genome Sequencing  - Bioinformatic II Course from Coursera
#
# Week 3 - Linear Spectrum problem
#    Input: An amino acid string Peptide
#    Output: The linear spectrum of Peptide
# Pseudocode:
"""
LinearSpectrum(Peptide, Alphabet, AminoAcidMass)
    PrefixMass(0) ← 0
    for i ← 1 to |Peptide|
        for every symbol s in Alphabet
            if s = i-th amino acid in Peptide
                PrefixMass(i) ← PrefixMass(i − 1) + AminoAcidMass[s]
    LinearSpectrum ← a list consisting of the single integer 0
    for i ← 0 to |Peptide| − 1
        for j ← i + 1 to |Peptide|
            add PrefixMass(j) − PrefixMass(i) to LinearSpectrum
    return sorted list LinearSpectrum
"""
#########################################################################################

# Constants
AMINO_ACID_MASS_MAP = {
    'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101,
    'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115, 'K': 128,
    'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156,
    'Y': 163, 'W': 186}

# linear Spectrum problem
def linear_spectrum(peptide_string):
    """
    Takes an peptide sequence to find its linear spectrum in ordered form

    Parameters:
    peptide_seq -> string of amino acid letters

    Returns:
    list of masses of subpeptides
    """

    if len(peptide_string) == 0:
        return [0]

    prefix_mass = [0]
    n = len(peptide_string)
    linear_spectrum = [0]

    for aa in peptide_string:
        mass = AMINO_ACID_MASS_MAP.get(aa)
        if mass is None:
            # Provide a clear error message to help debugging unexpected characters
            raise ValueError(f"Unknown amino acid '{aa}' in peptide '{peptide_string}'")
        prefix_mass.append(prefix_mass[-1] + mass)

    for i in range(n):
        for j in range(i+1, n+1):
            linear_spectrum.append(prefix_mass[j] - prefix_mass[i])

    return sorted(linear_spectrum)



###########################################################################

if __name__ == "__main__":

# Sample test
    peptide_string = "NQEL"
    # Expected answer = 0 113 114 128 129 242 242 257 370 371 484
    answer = linear_spectrum(peptide_string)
    # print(''.join(answer))
    print(*answer)


# # From file

#     # Get dataset
#     from pathlib import Path as partho
#     current_dir = partho(__file__).parent
#     filename = input("Please enter the filename: ")
#     file_path = current_dir / filename

#     #NOTE: Make sure to remove the extra line at the end of the file
#     with open(file_path, 'r') as file:
#        peptide_string = file.readline()

#     answer = linear_spectrum(peptide_string)
#     # print('\n'.join(answer))

#     with open("Wk3_4_output.txt", "w") as output_file:
#     #    output_file.write('\n'.join(answer))
#        output_file.write(" ".join(map(str, answer)))



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