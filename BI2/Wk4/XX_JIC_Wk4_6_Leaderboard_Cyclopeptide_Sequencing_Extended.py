#########################################################################################
# Genome Sequencing  - Bioinformatic II Course from Coursera
#
# Week 4 - Leaderboard Cyclopeptide Sequencing for Tyrocidine B1 - extended 57-200 AA
# Input: An integer N and a collection of integers Spectrum.
# Output: LeaderPeptide after running LeaderboardCyclopeptideSequencing(Spectrum, N).

"""
LeaderboardCyclopeptideSequencing(Spectrum, N)
    Leaderboard ← set containing only the empty peptide
    LeaderPeptide ← empty peptide
    while Leaderboard is non-empty
        Leaderboard ← Expand(Leaderboard)
        for each Peptide in Leaderboard
            if Mass(Peptide) = ParentMass(Spectrum)
                if Score(Peptide, Spectrum) > Score(LeaderPeptide, Spectrum)
                    LeaderPeptide ← Peptide
            else if Mass(Peptide) > ParentMass(Spectrum)
                remove Peptide from Leaderboard
        Leaderboard ← Trim(Leaderboard, Spectrum, N)
    output LeaderPeptide
"""
##### I am not going to import modules in this file as that is messing up the other files
##### as they use the full dict map, while here we just need the masses, not the aa letters
##### For this reason, many of the helper functions have been rewritten here
# NOTE: YOU HAVE TO RUN THIS FROM THE BI2 FOLDER WITH THE COMMAND:
# python -m Wk4.Wk4_4_Leaderboard_Cyclopeptide_Sequencing (NO .py)
#########################################################################################

# Imports:
from collections import Counter
from typing import List, Tuple

# Constants
AMINO_ACID_MASSES = [i for i in range(57, 201)]


######## Function 0 ########
# MAIN LEADERBOARD FUNCTION
#########################################

# The Parent Mass Logic: The parent mass (maximum value in experimental spectrum) represents the total mass of the
# complete Tyrocidine B1 molecule. Peptides below this are incomplete fragments still under construction.
# Only at parent_mass can you test if you've found the correct sequence by comparing cyclic spectra.
# The mass < parent_mass clause ensures growing peptides aren't discarded prematurely—they need more amino acids before evaluation.

def leaderboard_cyclopeptide_sequencing(experimental_spectrum: list[int], N: int):
    """
    Implements the Leaderboard Cyclopeptide Sequencing algorithm to find the best peptide

    Parameters:
        experimental_spectrum -> list of integers representing experimental spectrum
        N -> number of leaderboard peptides

    Returns:
        string of amino acid masses separated by '-'
    """

    # Use tuple-of-masses representation for peptides (mass-based) instead of amino-acid letters
    experimental_spectrum = sorted(experimental_spectrum)
    parent_mass = max(experimental_spectrum)
    leaderboard = [((), 0)]
    leader_peptides = []
    leader_score = -1

    while leaderboard:
        leaderboard = expand(leaderboard)

        # Dedup
        seen = set()
        new = []
        for peptide, mass in leaderboard:
            if peptide in seen:
                continue
            seen.add(peptide)
            new.append((peptide, mass))
        leaderboard = new
        # print("\n\nNew Leaderboard:", leaderboard)

        new_board = []
        for peptide, mass in leaderboard:
            # Critical: Parent mass indicates a complete peptide. Only full-length peptides are scored cyclically since they form rings in nature
            if mass == parent_mass:
                theo_cyclospectrum = cyclic_spectrum(peptide)
                score = spectrum_score(theo_cyclospectrum, experimental_spectrum)
                if score > leader_score:
                    leader_score = score
                    leader_peptides = [peptide]
                elif score == leader_score:
                    leader_peptides.append(peptide)
                new_board.append((peptide, mass))
            # Key insight: Incomplete peptides can't be scored cyclically (they're still linear chains). Keep them growing until they reach parent_mass
            elif mass < parent_mass:
                new_board.append((peptide, mass))
        leaderboard = trim(new_board, experimental_spectrum, N)
    return leader_peptides, leader_score


###########################################################################
############################# Helper Functions ############################
###########################################################################

######## Function 1 ########
# LINEAR SPECTRUM
# Adapted from Wk3.Wk3_4_LinearSpectrum
#########################################
def linear_spectrum(peptide: list[int]):
    """
    Takes an peptide sequence of masses to find its linear spectrum in ordered form

    Parameters:
        peptide -> list of amino acid masses (e.g. [113,128,186])

    Returns:
        list of masses of subpeptides
    """
    prefix_mass = [0]
    linear_spectrum = [0]
    n = len(peptide)

    if n == 0:
        return [0]

    for mass in peptide:
         prefix_mass.append(prefix_mass[-1] + mass)

    for i in range(n):
        for j in range(i+1, n+1):
            linear_spectrum.append(prefix_mass[j] - prefix_mass[i])

    return sorted(linear_spectrum)


######## Function 2 ########
# CYCLIC SPECTRUM
# Adapted from Wk3.Wk3_5_CyclicalSpectrum
#########################################
def cyclic_spectrum(peptide: list[int]) -> list:
    """
    Takes an peptide sequence to find its cyclic spectrum in ordered form

    Parameters:
        peptide -> list of amino acid masses (e.g. [113,128,186])

    Returns:
        list of masses of subpeptides
    """
    prefix_mass = [0]
    n = len(peptide)
    cyclic_spectrum = [0]

    for mass in peptide:
        prefix_mass.append(prefix_mass[-1] + mass)

    peptide_mass = prefix_mass[-1]
    for i in range(n):
        for j in range(i+1, n+1):
            cyclic_spectrum.append(prefix_mass[j] - prefix_mass[i])
            if i > 0 and j < n:
                cyclic_spectrum.append(peptide_mass - (prefix_mass[j] - prefix_mass[i]))
    return sorted(cyclic_spectrum)


####### Function 3 ########
# SPECTRUM SCORE
# Adapted from Wk4_2_LinearPeptide_Scoring
#########################################
def spectrum_score(theoretical_spectrum: list[int], experimental_spectrum: list[int]) -> int:
    """
    Compute the  score of a peptide represented as list of amino-acid masses (can be linear or cyclical)

    Parameters:
        theoretical_spectrum -> list of integer representing theoretical spectrum masses (e.g. [113,128,186])
        experimental_spectrum -> list of integers representing experimental spectrum

    Returns:
        integer score (counts multiplicities)
    """
    theo_count = Counter(theoretical_spectrum)
    exper_count = Counter(experimental_spectrum)
    score = 0
    for mass in theo_count:
        if mass in exper_count:
            score += min(theo_count[mass], exper_count[mass])
    return score


######## Function 4 ########
# EXPAND
# Adapted from Wk3.Wk3_8_CycloPeptide_Sequencing
##############################
def expand(peptides: list[int]) -> list:
    """Expand candidate peptides.

    Parameters:
        peptides: list of tuples of amino acid and mass e.g., [('X', 128), ]

    Returns:
        Expanded peptides with all other masses
    """
    unique_masses = set(AMINO_ACID_MASSES)
    expanded = []

    # Expand by appending integer masses
    for peptide, mass in peptides:
        for mass_val in unique_masses:
            expanded.append((peptide+(mass_val,), mass+mass_val))
            # print(expanded)
    return expanded


######## Function 5 ########
# Trim function
# Adapted from Wk4_3_Trim_Leaderboard
##############################

def trim(leaderboard: list[tuple[int, ...]], experimental_spectrum: list[int], N: int) -> list[tuple[int, ...]]:
    """
    Trims the leaderboard peptides (represented as tuples of masses) to the top N by linear score.
    Prevents exponential explosion of candidates. Poor-scoring partial peptides unlikely to become good complete ones

    Parameters:
        leaderboard -> list of peptides as tuples of masses
        experimental_spectrum -> list of integers representing experimental spectrum
        N -> number of leaderboard peptides to keep

    Returns:
        list of peptides (tuples) that remain after trimming (keeps ties at the N-th score)
    """
    # Compute (peptide, score) pairs
    scores = []
    for peptide, mass in leaderboard:
        theo_linear_spectrum = linear_spectrum(peptide)
        scores.append((spectrum_score(theo_linear_spectrum, experimental_spectrum), peptide, mass))

    # Sort by score descending
    scores.sort(key=lambda x: x[0], reverse=True)

    # If N is larger than the leaderboard size, keep all peptides sorted by score
    if N >= len(scores):
        return [(p, m) for s, p, m in scores]

    # Threshold score is the score of the N-th peptide (1-based). In 0-based index it's N-1
    threshold = scores[N-1][0]
    trimmed = [(p, m) for s, p, m in scores if s >= threshold]
    return trimmed


###########################################################################

if __name__ == "__main__":

# # Sample test
#     experimental_spectrum_raw = "0 97 99 114 128 147 147 163 186 227 241 242 244 260 261 262 283 291 333 340 357 385 389 390 390 405 430 430 447 485 487 503 504 518 543 544 552 575 577 584 632 650 651 671 672 690 691 738 745 747 770 778 779 804 818 819 820 835 837 875 892 917 932 932 933 934 965 982 989 1030 1039 1060  1061 1062 1078 1080 1081 1095 1136 1159 1175 1175 1194 1194 1208 1209 1223 1225 1322"
#     experimental_spectrum = list(map(int, experimental_spectrum_raw.split()))
#     N = 1000
#     # Expected answer = 34 linear peptides with max score of 87
#     leaders, best_score = leaderboard_cyclopeptide_sequencing(experimental_spectrum, N)
#     formatted = ["-".join(map(str, p)) for p in leaders]
#     print("Leader peptides are:\n")
#     print("\n".join(formatted))
#     print("\nBest score is:", best_score)
#     # print(*answer)

# From file

    # Get dataset
    from pathlib import Path as partho
    current_dir = partho(__file__).parent
    filename = input("Please enter the filename: ")
    file_path = current_dir / filename

    #NOTE: For Wk3_2, there is a newline character after the peptide string, so remove that
    with open(file_path, 'r') as file:
        # Read peptide string and remove ALL whitespace (spaces/newlines) — dataset may contain spaces
        # leaderboard = file.readline().split()
        N = int(file.readline())
        experimental_spectrum_raw = file.readline()
    experimental_spectrum = list(map(int, experimental_spectrum_raw.split()))

    leaders, best_score = leaderboard_cyclopeptide_sequencing(experimental_spectrum, N)
    formatted = ["-".join(map(str, p)) for p in leaders]
    answer = " ".join(formatted)

    with open("Wk4_6_output.txt", "w") as output_file:
        # output_file.write(str(answer))
        output_file.write(answer)


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