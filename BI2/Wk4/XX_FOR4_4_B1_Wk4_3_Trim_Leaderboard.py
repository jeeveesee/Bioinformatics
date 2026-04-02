#########################################################################################
# Genome Sequencing  - Bioinformatic II Course from Coursera
#
# Week 4 - Trimming algorithm for leaderboard scoring
# Input: A collection of peptides Leaderboard, a collection of integers Spectrum, and an integer N.
# Output: The N highest-scoring linear peptides on Leaderboard with respect to Spectrum.
# Step 2: Create the trimming logic for the leaderboard
#   - Find linear score
#   - Order them in descending order
#   - Find the Top N linear scores along with multiplicities at the end and then voila!
# This is part 2 of the Trim Leaderboard algorithm
# ONLY the Trim function uses linear scoring. The actual cyclopeptide sequencing uses cycloscoring

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
# NOTE: YOU HAVE TO RUN THIS FROM THE BI2 FOLDER WITH THE COMMAND:
# python -m Wk4.Wk4_3_Trim_Leaderboard (NO .py)
#########################################################################################

# Imports:
from Wk4.Wk4_2_LinearPeptide_Scoring import linear_scoring


def trim(leaderboard: list[tuple[int, ...]], experimental_spectrum: list[int], N: int) -> list[tuple[int, ...]]:
    """
    Trims the leaderboard peptides (represented as tuples of masses) to the top N by linear score.

    Parameters:
    leaderboard -> list of peptides as tuples of masses
    experimental_spectrum -> list of integers representing experimental spectrum
    N -> number of leaderboard peptides to keep

    Returns:
    list of peptides (tuples) that remain after trimming (keeps ties at the N-th score)
    """

    # Compute (peptide, score) pairs
    scores = []
    for peptide in leaderboard:
        scores.append((peptide, linear_scoring(peptide, experimental_spectrum)))

    # Sort by score descending
    scores.sort(key=lambda x: x[1], reverse=True)

    # If N is larger than the leaderboard size, keep all peptides sorted by score
    if N >= len(scores):
        return [p for p, s in scores]

    # Threshold score is the score of the N-th peptide (1-based). In 0-based index it's N-1
    threshold = scores[N-1][1]

    trimmed = [p for p, s in scores if s >= threshold]
    return trimmed

    # The following works but too confusing
    ##########
    # # Find the first index after the top-N where score < threshold (stop early for efficiency)
    # cutoff = len(scores)
    # for i in range(N, len(scores)):
    #     if scores[i][1] < threshold:
    #         cutoff = i
    #         break

    # # Return peptides up to cutoff (this keeps ties at the N-th position)
    # return [p for p, s in scores[:cutoff]]
    ##########

    ##########
    # The following should also work as per the text but is kinda annoying and too verbose
    # for j in range(N, len_leaderboard):
    #     if scores[j] < scores[N]:
    #         leaderboard = leaderboard[:j-1]
    #         return leaderboard
    # return leaderboard
    ##########



###########################################################################

if __name__ == "__main__":

# Sample test
    leaderboard = ['LAST', 'ALST', 'TLLT', 'TQAS']
    experimental_spectrum = [0, 71, 87, 101, 113, 158, 184, 188, 259, 271, 372]
    N = 2
    # Expected answer = LAST ALST
    answer = trim(leaderboard, experimental_spectrum, N)
    # print('\n'.join(answer))
    print(*answer)

# # From file

#     # Get dataset
#     from pathlib import Path as partho
#     current_dir = partho(__file__).parent
#     filename = input("Please enter the filename: ")
#     file_path = current_dir / filename

#     #NOTE: For Wk3_2, there is a newline character after the peptide string, so remove that
#     with open(file_path, 'r') as file:
#         # Read peptide string and remove ALL whitespace (spaces/newlines) — dataset may contain spaces
#         leaderboard = file.readline().split()
#         experimental_spectrum_raw = file.readline()
#         N = int(file.readline())

#     experimental_spectrum = list(map(int, experimental_spectrum_raw.split()))

#     answer = trim(leaderboard, experimental_spectrum, N)
#     # print(' '.join(answer))
#     print(*answer)

#     with open("Wk4/Wk4_3_output.txt", "w") as output_file:
#         # output_file.write(str(answer))
#         output_file.write(' '.join(answer))


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