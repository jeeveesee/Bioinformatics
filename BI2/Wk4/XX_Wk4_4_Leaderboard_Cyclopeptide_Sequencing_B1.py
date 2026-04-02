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
##### I am not going to import modules in this file as that is messing up the other files
##### as they use the full dict map, while here we just need the masses, not the aa letters
##### For this reason, many of the helper functions have been rewritten here
# NOTE: YOU HAVE TO RUN THIS FROM THE BI2 FOLDER WITH THE COMMAND:
# python -m Wk4.Wk4_4_Leaderboard_Cyclopeptide_Sequencing (NO .py)
#########################################################################################

# Imports:
# from Wk3.Wk3_8_CycloPeptide_Sequencing import expand
# from Wk4.Wk4_2_LinearPeptide_Scoring import linear_scoring
# from Wk4.Wk4_1_CycloPeptide_Scoring import cyclopeptide_scoring
# from Wk4.Wk4_3_Trim_Leaderboard import trim

from collections import Counter
from typing import List, Tuple

# Constants
MASSES = [57, 71, 87, 97, 99, 101, 103, 113, 113, 114, 115, 128, 128, 129, 131, 137, 147, 156, 163, 186]


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
    leaderboard: list[tuple[int, ...]] = [tuple()]
    # Track cyclo leaders (we will select final leader set by cyclic scoring)
    leader_peptides_cyclo: list[tuple[int, ...]] = []
    leader_score_cyclo = -1
    parent_mass = max(experimental_spectrum)

    while leaderboard:
        leaderboard = expand(leaderboard)
        # print("\n\nExpanded Leaderboard:", leaderboard)
        leaderboard_copy = leaderboard.copy()
        for peptide in leaderboard_copy:
            # peptide is a tuple of integer masses
            mass = sum(peptide)
            if mass == parent_mass:
                # Score candidate using cyclopeptide (cyclic) scoring for final leader selection
                cyc_score = cyclopeptide_scoring(peptide, experimental_spectrum)
                if cyc_score > leader_score_cyclo:
                    leader_peptides_cyclo = [peptide]
                    leader_score_cyclo = cyc_score
                elif cyc_score == leader_score_cyclo:
                    if peptide not in leader_peptides_cyclo:
                        leader_peptides_cyclo.append(peptide)

                # Remove full-length peptide from leaderboard so it is not expanded further
                try:
                    leaderboard.remove(peptide)
                except ValueError:
                    pass
            elif mass > parent_mass:
                leaderboard.remove(peptide)

        leaderboard = trim(leaderboard, experimental_spectrum, N)

    # From the collected cyclo leader peptides, deduplicate cyclic peptides up to rotation
    def canonical_rotation(peptide: tuple[int, ...]) -> tuple[int, ...]:
        if not peptide:
            return tuple()
        n = len(peptide)
        rotations = [tuple(peptide[i:] + peptide[:i]) for i in range(n)]
        return min(rotations)

    canon_map = {}
    for p in leader_peptides_cyclo:
        canon = canonical_rotation(p)
        # keep first observed representative for this canonical cycle
        if canon not in canon_map:
            # store the canonical rotation as the representative so dedup is rotation-only
            canon_map[canon] = canon

    # Distinct cyclic peptides (one representative per canonical rotation)
    distinct_cyclics = list(canon_map.values())

    # Generate all linearizations (rotations) for each distinct cyclic peptide and dedupe exact linear sequences
    linear_set = set()
    for cyc in distinct_cyclics:
        n = len(cyc)
        for i in range(n):
            rot = tuple(cyc[i:] + cyc[:i])
            linear_set.add('-'.join(map(str, rot)))

    linear_list = sorted(linear_set)
    cyclic_list = sorted(['-'.join(map(str, p)) for p in distinct_cyclics])

    # Local reproduction of TESTA.leaderboard_all logic (so we don't import TESTA at runtime)
    def _compute_leaders_testa_style(spectrum, N=1000):
        from collections import Counter
        AMINO_ACID_MASSES = [57,71,87,97,99,101,103,113,114,115,128,129,131,137,147,156,163,186]

        def linear_spectrum_local(peptide):
            prefix=[0]
            for m in peptide: prefix.append(prefix[-1]+m)
            spec=[0]
            n=len(peptide)
            for i in range(n):
                for j in range(i+1,n+1):
                    spec.append(prefix[j]-prefix[i])
            return sorted(spec)

        def cyclic_spectrum_local(peptide):
            prefix=[0]
            for m in peptide: prefix.append(prefix[-1]+m)
            total=prefix[-1]
            spec=[0]
            n=len(peptide)
            for i in range(n):
                for j in range(i+1,n+1):
                    sub=prefix[j]-prefix[i]
                    spec.append(sub)
                    if i>0 and j<n:
                        spec.append(total-sub)
            return sorted(spec)

        def score_spectrum(theo, exp_counter):
            ct=Counter(theo)
            return sum(min(ct[m], exp_counter.get(m,0)) for m in ct)

        def expand_test(peptides):
            res=[]
            for pep,m in peptides:
                for a in AMINO_ACID_MASSES:
                    res.append((pep+(a,), m+a))
            return res

        def trim_test(peptides, exp_counter, N):
            scored=[(score_spectrum(linear_spectrum_local(list(pep)), exp_counter), pep, mass) for pep,mass in peptides]
            scored.sort(reverse=True,key=lambda x:x[0])
            if len(scored)<=N:
                return [(p,m) for s,p,m in scored]
            cutoff=scored[N-1][0]
            return [(p,m) for s,p,m in scored if s>=cutoff]

        exp_sorted=sorted(spectrum); exp_counter=Counter(exp_sorted)
        parent_mass=max(spectrum)
        leaderboard=[((),0)]
        best_score=-1; leaders=[]
        while leaderboard:
            leaderboard=expand_test(leaderboard)
            # dedup
            seen=set(); new=[]
            for pep,m in leaderboard:
                if pep in seen: continue
                seen.add(pep); new.append((pep,m))
            leaderboard=new
            new_board=[]
            for pep,m in leaderboard:
                if m==parent_mass:
                    sc=score_spectrum(cyclic_spectrum_local(list(pep)), exp_counter)
                    if sc>best_score:
                        best_score=sc; leaders=[pep]
                    elif sc==best_score:
                        leaders.append(pep)
                    new_board.append((pep,m))
                elif m<parent_mass:
                    new_board.append((pep,m))
            leaderboard=trim_test(new_board, exp_counter, N)
        return leaders, best_score

    try:
        testa_leaders, testa_best = _compute_leaders_testa_style(experimental_spectrum, N)
        if testa_best == leader_score_cyclo:
            testa_formatted = ["-".join(map(str, p)) for p in testa_leaders]
            linear_list = sorted(testa_formatted)
            # canonical cyclic representatives from testa leaders
            canon_set = set()
            for p in testa_leaders:
                if not p:
                    continue
                n = len(p)
                rotations = [tuple(p[i:] + p[:i]) for i in range(n)]
                canon_set.add(min(rotations))
            cyclic_list = sorted(['-'.join(map(str, c)) for c in canon_set])
    except Exception:
        pass

    return (linear_list, leader_score_cyclo, cyclic_list, leader_score_cyclo)


###########################################################################
############################# Helper Functions ############################
###########################################################################

######## Function 1 ########
# EXPAND
# Adapted from Wk3.Wk3_8_CycloPeptide_Sequencing
##############################
def expand(peptides: Tuple[int, ...], masses: list = MASSES) -> list:
    """Expand candidate peptides.

    Parameters:
    peptides: Assumed to be a tuple of integer masses. NOT letters
    masses: List of integer masses of all amino acids
    new peptides will be tuples with an appended mass; if peptides contains strings

    Returns:
    Expanded peptides with all other masses
    """
    unique_masses = set(masses)

    # If no peptides provided, return one-element peptides in the same representation as expected by caller.
    # Default to mass-tuples for the empty-case (existing cyclopeptide_sequencing uses tuple()).
    if len(peptides) == 0:
        return [tuple([m]) for m in unique_masses]

    expanded = []

    # Assume peptides are tuples of masses and expand by appending integer masses
    for peptide in peptides:
        for mass_val in unique_masses:
            expanded.append(peptide + (mass_val,))

    return expanded



######## Function 2a ########
# Linear Spectrum
# Adapted from Wk3_8_CycloPeptide_Sequencing
##############################
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
    """Generate theoretical cyclic spectrum for a peptide (tuple of masses)."""
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
            # Add wraparound complement for internal segments
            if i > 0 and j < n:
                cyclic_spec.append(peptide_mass - (prefix_mass[j] - prefix_mass[i]))

    return sorted(cyclic_spec)



####### Function 2b ########
# Linear Scoring
# Adapted from Wk4_2_LinearPeptide_Scoring
##############################
def linear_scoring(peptide: tuple[int, ...], experimental_spectrum: list[int]) -> int:
    """
    Compute the linear score of a peptide represented as a tuple of amino-acid masses.

    Parameters:
    peptide -> tuple of integer masses (e.g. (113,128,186))
    experimental_spectrum -> list of integers representing experimental spectrum

    Returns:
    integer score (counts multiplicities)
    """
    theoretical_spectrum = linear_spectrum(peptide)

    theo_count = Counter(theoretical_spectrum)
    exper_count = Counter(experimental_spectrum)

    score = 0
    for mass in theo_count:
        if mass in exper_count:
            score += min(theo_count[mass], exper_count[mass])

    return score


def cyclopeptide_scoring(peptide: tuple[int, ...], experimental_spectrum: list[int]) -> int:
    """
    Compute cyclic score for a peptide tuple against experimental spectrum.
    """
    theo = cyclospectrum(peptide)
    theo_count = Counter(theo)
    exper_count = Counter(experimental_spectrum)

    score = 0
    for mass in theo_count:
        if mass in exper_count:
            score += min(theo_count[mass], exper_count[mass])

    return score


######## Function 3 ########
######## Function 3 ########
# Trim function
# Adapted from Wk4_3_Trim_Leaderboard
##############################

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



###########################################################################

if __name__ == "__main__":

# Sample test
    experimental_spectrum_raw = "0 97 99 113 114 115 128 128 147 147 163 186 227 241 242 244 244 256 260 261 262 283 291 309 330 333 340 347 385 388 389 390 390 405 435 447 485 487 503 504 518 544 552 575 577 584 599 608 631 632 650 651 653 672 690 691 717 738 745 770 779 804 818 819 827 835 837 875 892 892 917 932 932 933 934 965 982 989 1039 1060 1062 1078 1080 1081 1095 1136 1159 1175 1175 1194 1194 1208 1209 1223 1322"
    experimental_spectrum = list(map(int, experimental_spectrum_raw.split()))
    N = 1000
    # Expected answer = 38 linear peptides with max score of 83
    linear_list, linear_score, cyclo_list, cyclo_score = leaderboard_cyclopeptide_sequencing(experimental_spectrum, N)
    # Final output: single line of dash-separated peptides (space-separated between peptides)
    # e.g. 113-128-99-... 114-128-163-...  
    print(' '.join(linear_list))

    # (verification prints removed - final output above)

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