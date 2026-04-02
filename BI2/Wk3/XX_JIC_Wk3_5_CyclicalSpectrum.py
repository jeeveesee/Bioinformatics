#########################################################################################
# Genome Sequencing  - Bioinformatic II Course from Coursera
#
# Week 3 - Cyclic Spectrum problem
#    Input: An amino acid string Peptide
#    Output: The cyclic spectrum of Peptide
# Pseudocode:
"""
CyclicSpectrum(Peptide, Alphabet, AminoAcidMass)
    PrefixMass(0) ← 0
    for i ← 1 to |Peptide|
        for every symbol s in Alphabet
            if s = i-th amino acid in Peptide
                PrefixMass(i) ← PrefixMass(i - 1) + AminoAcidMass[s]
    peptideMass ← PrefixMass(|Peptide|)
    CyclicSpectrum ← a list consisting of the single integer 0
    for i ← 0 to |Peptide| - 1
        for j ← i + 1 to |Peptide|
            add PrefixMass(j) - PrefixMass(i) to CyclicSpectrum
            if i > 0 and j < |Peptide|
                add peptideMass - (PrefixMass(j) - PrefixMass(i)) to CyclicSpectrum
    return sorted list CyclicSpectrum
"""
#########################################################################################

# Constants
AMINO_ACID_MASS_MAP = {
    'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101,
    'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115, 'K': 128,
    'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156,
    'Y': 163, 'W': 186}

# linear Spectrum problem
def cyclic_spectrum(amino_acid_mass_map, peptide_string):
    """
    Takes an peptide sequence to find its cyclic spectrum in ordered form

    Parameters:
    peptide_seq -> string of amino acid letters
    amino_acid_mass_map -> dict of mass of each amino acid letter

    Returns:
    list of masses of subpeptides
    """
    prefix_mass = [0]
    n = len(peptide_string)
    cyclic_spectrum = [0]

    for idx, aa in enumerate(peptide_string):
        mass = amino_acid_mass_map.get(aa)
        if mass is None:
            raise ValueError(f"Unknown amino acid '{aa}' at position {idx} in peptide_string (maybe contains newline or whitespace). Peptide: {repr(peptide_string)}")
        prefix_mass.append(prefix_mass[-1] + mass)

    peptide_mass = prefix_mass[-1]

    for i in range(n):
        for j in range(i+1, n+1):
            cyclic_spectrum.append(prefix_mass[j] - prefix_mass[i])
            if i > 0 and j < n:
                cyclic_spectrum.append(peptide_mass - (prefix_mass[j] - prefix_mass[i]))

    return sorted(cyclic_spectrum)



###########################################################################

if __name__ == "__main__":

# Sample test
    peptide_string = "RESINFNRGLVFYNRANWEPQSRGRDCTSRHMHNAWY" #"LEQN"
    # Expected answer = 0 113 114 128 129 227 242 242 257 355 356 370 371 484
    answer = cyclic_spectrum(AMINO_ACID_MASS_MAP, peptide_string)
    # print(''.join(answer))
    print(*answer)
    # For exam
    # print(' '.join(map(str,answer)) in " 0 71 101 113 131 184 202 214 232 285 303 315 345 416")


# # # From file

#     # Get dataset
#     from pathlib import Path as partho
#     current_dir = partho(__file__).parent
#     filename = input("Please enter the filename: ")
#     file_path = current_dir / filename

#     #NOTE: Make sure to remove the extra line at the end of the file
#     with open(file_path, 'r') as file:
#        peptide_string = file.readline()

#     answer = cyclic_spectrum(AMINO_ACID_MASS_MAP, peptide_string)
#     # print('\n'.join(answer))

#     with open("Wk3_5_output.txt", "w") as output_file:
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