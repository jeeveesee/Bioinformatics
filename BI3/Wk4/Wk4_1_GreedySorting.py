#########################################################################################
# Genome Sequencing  - Bioinformatic III Course from Coursera
#
# Week 4 - Random Breakage Model of Chromosome Evolution - Greedy Sorting
# Code Challenge: Implement GreedySorting.
# Input: A permutation P.
# Output: The sequence of permutations corresponding to applying GreedySorting to P,
# ending with the identity permutation.
#########################################################################################


def format_permutation(perm):
    """
    Formats a permutation (list of signed integers) as a space-separated string.

    Parameters:
    perm -> list of signed integers e.g. [-1, -4, 3, 5, -2]

    Returns:
    space-separated string with explicit + or - signs e.g. "-1 -4 +3 +5 -2"
    """
    result = []
    for x in perm:
        if x > 0:
            result.append(f"+{x}")
        else:
            result.append(str(x))
    return " ".join(result)


def reverse_segment(perm, i, j):
    """
    Reverses the segment of perm from index i to j (inclusive) and flips the sign
    of every element in that segment — as required by signed permutation reversals.

    Parameters:
    perm -> list of signed integers
    i    -> start index of the reversal (inclusive)
    j    -> end index of the reversal (inclusive)

    Returns:
    new permutation with the segment reversed and all signs in that segment negated
    """
    return perm[:i] + [-x for x in reversed(perm[i : j + 1])] + perm[j + 1 :]


def greedy_sorting(permutation):
    """
    Implements GreedySorting to sort a signed permutation into the identity permutation.
    At each position k it finds element k, brings it to position k (at most 1 reversal),
    then fixes the sign if needed (at most 1 more reversal). Outputs every intermediate step.

    Parameters:
    permutation -> list of strings with +/- signs e.g. ["-3", "+4", "+1", "+5", "-2"]

    Returns:
    list of strings, each being a space-separated permutation produced during sorting
    """
    # Convert string elements like "+4" / "-3" to signed integers
    perm = [int(x) for x in permutation]
    n = len(perm)
    steps = []

    for k in range(n):  # for k ← 1 to |P|  (0-indexed, so element = k+1)
        element = k + 1

        if perm[k] == element:  # already in place with correct sign — skip
            continue

        # Find the position of +element or -element, starting from index k
        for i in range(k, n):
            if abs(perm[i]) == element:
                pos = i
                break

        # Step 1: if element is not at index k, reverse the segment to bring it there
        if pos != k:
            perm = reverse_segment(perm, k, pos)
            steps.append(format_permutation(perm))

        # Step 2: if element arrived with the wrong sign, flip it (single-element reversal)
        if perm[k] == -element:
            perm = reverse_segment(perm, k, k)
            steps.append(format_permutation(perm))

    return steps


###########################################################################

if __name__ == "__main__":
    # # Sample test
    # permutation = [-3, 4, 1, 5, -2]
    # # Expected answer =
    # # -1 -4 +3 +5 -2
    # # +1 -4 +3 +5 -2
    # # +1 +2 -5 -3 +4
    # # +1 +2 +3 +5 +4
    # # +1 +2 +3 -4 -5
    # # +1 +2 +3 +4 -5
    # # +1 +2 +3 +4 +5
    # answer = greedy_sorting(permutation)
    # print("\n".join(answer))
    # # print(answer)

    # # # From file

    # # Get dataset
    # from pathlib import Path as partho

    # current_dir = partho(__file__).parent
    # filename = input("Please enter the filename: ")
    # file_path = current_dir / filename

    # # NOTE: For Wk3_2, there is a newline character after the peptide string, so remove that
    # with open(file_path, "r") as file:
    #     permutation = list(map(int, file.read().split()))

    # # print(f"{permutation=}")
    # answer = greedy_sorting(permutation)
    # # print("\n".join(answer))
    # # print(answer)

    # with open("Wk4_1_output.txt", "w") as output_file:
    #     output_file.write("\n".join(answer))

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

    # For exam
    permutation = [20, 7, 10, 9, 11, 13, 18, -8, -6,
                   -14, 2, -4, -16, 15, 1, 17, 12, -5, 3, -19]
    # # Expected answer =
    # # -1 -4 +3 +5 -2
    # # +1 -4 +3 +5 -2
    # # +1 +2 -5 -3 +4
    # # +1 +2 +3 +5 +4
    # # +1 +2 +3 -4 -5
    # # +1 +2 +3 +4 -5
    # # +1 +2 +3 +4 +5
    answer = greedy_sorting(permutation)
    print("\n".join(answer))
    # print(answer)
