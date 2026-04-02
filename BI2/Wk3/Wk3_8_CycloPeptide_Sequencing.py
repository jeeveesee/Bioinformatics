#########################################################################################
# Genome Sequencing  - Bioinformatic II Course from Coursera
#
# Week 3 - CycloPeptide Sequencing
# Input: String of integers
# Output: String of peptides

"""
CyclopeptideSequencing(Spectrum)
    CandidatePeptides ← a set containing only the empty peptide FinalPeptides ← empty list of strings
    while CandidatePeptides is nonempty
        CandidatePeptides ← Expand(CandidatePeptides)
        for each peptide Peptide in CandidatePeptides
            if Mass(Peptide) = ParentMass(Spectrum)
                if Cyclospectrum(Peptide) = Spectrum and Peptide is not in FinalPeptides
                    append Peptide to FinalPeptides
                remove Peptide from CandidatePeptides
            else if Peptide is not consistent with Spectrum
                remove Peptide from CandidatePeptides
    return FinalPeptides
"""

# Linear_Spectrum: Couldn't use Wk3_4 Linear Spectrum function because that is based on
#                  the peptide letter e.g., PVKF, but not really the mass.
#                  So had to rewrite. But, it's pretty much the same

# From Phillip Compeau:
# Remember: when checking consistency, check consistency as a linear peptide. 
# When checking if the spectrum of a peptide is equal to the given spectrum, use the cyclic spectrum of the peptide.
#########################################################################################

# Imports
from collections import Counter
from typing import List, Tuple

# Constants
AMINO_ACID_MASS_MAP = {
    'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101,
    'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115, 'K': 128,
    'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156,
    'Y': 163, 'W': 186}

########################################################################################
########################################################################################
# Helper functiions
########################################################################################
########################################################################################

def mass(peptide: Tuple[int, ...]) -> int:
    """Calculate the mass of a peptide (tuple of amino acid masses)"""
    return sum(peptide)

def parent_mass(spectrum: List[int]) -> int:
    """Get parent mass (largest value in spectrum)"""
    return max(spectrum)

def linear_spectrum(peptide: Tuple[int, ...]) -> List[int]:
    """Generate theoretical linear spectrum"""
    if len(peptide) == 0:
        return [0]

    prefix_mass = [0]
    n = len(peptide)
    linear_spectrum = [0]

    for mass_val in peptide:
        prefix_mass.append(prefix_mass[-1] + mass_val)

    for i in range(n):
        for j in range(i+1, n+1):
            linear_spectrum.append(prefix_mass[j] - prefix_mass[i])

    return sorted(linear_spectrum)

def cyclospectrum(peptide: Tuple[int, ...]) -> List[int]:
    """Generate theoretical cyclic spectrum"""

    prefix_mass = [0]
    n = len(peptide)

    if n == 0:
        return [0]

    for mass_val in peptide:
        prefix_mass.append(prefix_mass[-1] + mass_val)

    peptide_mass = prefix_mass[-1]
    cyclic_spec = [0]

    for i in range(n):
        for j in range(i + 1, n + 1):
            cyclic_spec.append(prefix_mass[j] - prefix_mass[i])
            # Add wraparound subpeptide (complement) for internal segments
            if i > 0 and j < n:
                cyclic_spec.append(peptide_mass - (prefix_mass[j] - prefix_mass[i]))

    return sorted(cyclic_spec)

def is_consistent(peptide: Tuple[int, ...], spectrum: List[int]) -> bool:
    """Check if peptide is consistent with spectrum"""
    linear_spec = linear_spectrum(peptide)

    # Create frequency counters
    spectrum_freq = Counter(spectrum)
    linear_freq = Counter(linear_spec)

    # Check if all masses in linear spectrum are in experimental spectrum
    # with sufficient frequency
    for mass_val, count in linear_freq.items():
        if spectrum_freq[mass_val] < count:
            return False

    return True

def expand(peptides: list, masses: dict = AMINO_ACID_MASS_MAP) -> list:
    """Expand candidate peptides.

    This function supports two peptide representations:
    - tuple of integer masses (used by the mass-based cyclopeptide functions in this file)
    - string of amino-acid letters (used by other modules that operate on letter-based peptides)

    The expansion preserves the input representation: if peptides contains tuples,
    new peptides will be tuples with an appended mass; if peptides contains strings,
    new peptides will be strings with an appended amino-acid letter.
    """
    unique_masses = set(masses.values())

    # If no peptides provided, return one-element peptides in the same representation as expected by caller.
    # Default to mass-tuples for the empty-case (existing cyclopeptide_sequencing uses tuple()).
    if len(peptides) == 0:
        return [tuple([m]) for m in unique_masses]

    expanded = []

    sample = peptides[0]
    # If peptides are strings (amino-acid letters), expand by appending letters
    if isinstance(sample, str):
        for peptide in peptides:
            for aa in masses.keys():
                expanded.append(peptide + aa)
        return expanded

    # Otherwise assume peptides are tuples of masses and expand by appending integer masses
    for peptide in peptides:
        for mass_val in unique_masses:
            expanded.append(peptide + (mass_val,))

    return expanded

def peptide_to_string(peptide: Tuple[int, ...]) -> str:
    """Convert peptide (mass tuple) to string representation"""
    return '-'.join(map(str, peptide))

def cyclopeptide_sequencing(spectrum: List[int]) -> List[str]:
    """
    Main cyclopeptide sequencing algorithm
    Returns list of peptide strings that match the spectrum
    """
    candidate_peptides = [tuple()]  # Start with empty peptide
    final_peptides = []
    parent = parent_mass(spectrum)

    while candidate_peptides:
        # Expand all candidate peptides
        candidate_peptides = expand(candidate_peptides)

        # Filter peptides
        new_candidates = []
        for peptide in candidate_peptides:
            peptide_mass = mass(peptide)

            if peptide_mass == parent:
                # Check if cyclospectrum matches
                if cyclospectrum(peptide) == spectrum:
                    peptide_str = peptide_to_string(peptide)
                    if peptide_str not in final_peptides:
                        final_peptides.append(peptide_str)
                # Don't add to new_candidates (remove from search)
            elif is_consistent(peptide, spectrum):
                # Keep peptide if it's consistent
                new_candidates.append(peptide)
            # else: remove from search (not consistent)

        candidate_peptides = new_candidates

    return final_peptides


###########################################################################

if __name__ == "__main__":

# Sample test
    spectrum_raw = "0 113 128 186 241 299 314 427"
    # spectrum_raw = "0 71 101 113 131 184 202 214 232 285 303 315 345 416" # from exam
    spectrum = [int(num) for num in spectrum_raw.split()]
    # Expected answer = "186-128-113 186-113-128 128-186-113 128-113-186 113-186-128 113-128-186"
    answer = cyclopeptide_sequencing(spectrum)
    ## print(''.join(answer))
    print(*answer)


# # # From file

#     # Get dataset
#     from pathlib import Path as partho
#     current_dir = partho(__file__).parent
#     filename = input("Please enter the filename: ")
#     file_path = current_dir / filename

#     #NOTE: Make sure to remove the extra line at the end of the file
#     with open(file_path, 'r') as file:
#        spectrum_raw = file.readline()
#     spectrum = [int(num) for num in spectrum_raw.split()]

#     answer = cyclopeptide_sequencing(spectrum)
#     # print('\n'.join(answer))

#     with open("Wk3_8_output.txt", "w") as output_file:
#     #   output_file.write('\n'.join(answer))
#       output_file.write(" ".join(map(str, answer)))
#     #   output_file.write(str(answer))



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