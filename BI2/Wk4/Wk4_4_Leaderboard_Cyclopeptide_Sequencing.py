#########################################################################################
# Genome Sequencing  - Bioinformatic II Course from Coursera
#
# Week 4 - Leaderboard Cyclopeptide Sequencing
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
# Btw, mass should be equal to parent mass because only those peptides
# can be candidates for the spectrum. Others will be too light/partial
# NOTE: YOU HAVE TO RUN THIS FROM THE BI2 FOLDER WITH THE COMMAND:
# python -m Wk4.Wk4_4_Leaderboard_Cyclopeptide_Sequencing (NO .py)
#########################################################################################

# Imports:
from Wk3.Wk3_8_CycloPeptide_Sequencing import expand
from Wk4.Wk4_2_LinearPeptide_Scoring import linear_scoring
from Wk4.Wk4_1_CycloPeptide_Scoring import cyclopeptide_scoring
from Wk4.Wk4_3_Trim_Leaderboard import trim

# Constants
AMINO_ACID_MASS_MAP = {
    'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101,
    'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115, 'K': 128,
    'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156,
    'Y': 163, 'W': 186}

def leaderboard_cyclopeptide_sequencing(experimental_spectrum: list[int], N: int):
    """
    Implements the Leaderboard Cyclopeptide Sequencing algorithm to find the best peptide

    Parameters:
    experimental_spectrum -> list of integers representing experimental spectrum
    N -> number of leaderboard peptides

    Returns:
    string of amino acid masses separated by '-'
    """

    leaderboard = ['']
    leader_peptide = ''
    parent_mass = max(experimental_spectrum)

    while leaderboard:
        leaderboard = expand(leaderboard)
        # print("\n\nExpanded Leaderboard:", leaderboard)
        leaderboard_copy = leaderboard.copy()
        for peptide in leaderboard_copy:
            mass = sum(AMINO_ACID_MASS_MAP[aa] for aa in peptide)
            if mass == parent_mass:
                if cyclopeptide_scoring(peptide, experimental_spectrum) > cyclopeptide_scoring(leader_peptide, experimental_spectrum):
                    leader_peptide = peptide
            elif mass > parent_mass:
                leaderboard.remove(peptide)
        leaderboard = trim(leaderboard, experimental_spectrum, N)

    # Convert leader_peptide to mass string
    leader_peptide_masses = [str(AMINO_ACID_MASS_MAP[aa]) for aa in leader_peptide]
    return leader_peptide_masses



###########################################################################

if __name__ == "__main__":

# Sample test
    experimental_spectrum = [0, 71, 113, 129, 147, 200, 218, 260, 313, 331, 347, 389, 460]
    N = 10
    # Expected answer = 113-147-71-129
    # There might be multiple solutions
    answer = leaderboard_cyclopeptide_sequencing(experimental_spectrum, N)
    print('-'.join(answer))
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
#         N = int(file.readline())
#         experimental_spectrum_raw = file.readline()

#     experimental_spectrum = list(map(int, experimental_spectrum_raw.split()))

#     answer = leaderboard_cyclopeptide_sequencing(experimental_spectrum, N)
#     print('-'.join(answer))
#     # print(*answer)

#     with open("Wk4/Wk4_4_output.txt", "w") as output_file:
#         # output_file.write(str(answer))
#         output_file.write('-'.join(answer))


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