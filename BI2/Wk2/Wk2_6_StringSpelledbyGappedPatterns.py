#########################################################################################
# Genome Sequencing  - Bioinformatic II Course from Coursera
#
# Week 2 - Strings spelled by gapped patterns
# This algo generalizes generalizes the approach to an arbitrary sequence GappedPatterns of (k, d)-mers.
# It constructs strings PrefixString and SuffixString as described above,
# and checks whether they have perfect overlap (i.e., form the prefix and suffix of a
# reconstructed string). It also assumes that the number of (k, d)-mers in GappedPatterns
# is at least d; otherwise, it is impossible to reconstruct a contiguous string.

# pseudocode:
"""
StringSpelledByGappedPatterns(GappedPatterns, k, d)
    FirstPatterns ← the sequence of initial k-mers from GappedPatterns
    SecondPatterns ← the sequence of terminal k-mers from GappedPatterns
    PrefixString ← StringSpelledByPatterns(FirstPatterns, k)
    SuffixString ← StringSpelledByPatterns(SecondPatterns, k)
    for i = k + d + 1 to |PrefixString|
        if the i-th symbol in PrefixString does not equal the (i - k - d)-th symbol in SuffixString
            return "there is no string spelled by the gapped patterns"
    return PrefixString concatenated with the last k + d symbols of SuffixString
"""
# Input:
# Integers k and d followed by a sequence of (k, d)-mers (a1|b1), … , (an|bn)
# such that Suffix(ai|bi) = Prefix(ai+1|bi+1) for 1 ≤ i ≤ n-1.
# Output:
# A string Text of length k + d + k + n - 1 such that the i-th (k, d)-mer in Text
# is equal to (ai|bi)  for 1 ≤ i ≤ n (if such a string exists).
#########################################################################################
#########################################################################################

# Don't need this - I did this for myself to learn
def string_spelled_by_patterns(patterns, k):
    """
    Recreates full genome from a string of patterns,
    where all the inputs are sanitized to be in the right format
    :param patterns -> patterns as a string separated by a space
    :param k -> length of kmers
    :return -> single genome
    """
    # Clean patterns for right format
    patterns_sanitized = patterns.split()

    # Create genome from patterns
    genome = patterns_sanitized[0]
    for pattern in patterns_sanitized[1:]:
        genome += pattern[-1]

    return genome

######################
# Actual model
######################
def string_spelled_by_gapped_patterns(gapped_patterns, k, d):
    """
    Recreates full genome from a string of gapped patterns,
    where all the inputs are sanitized to be in the right format
    :param gapped_patterns -> gapped patterns as a string separated by a space
    :param k -> length of kmers
    :param d -> gap between kmers
    :return -> single genome
    """
    # Clean gapped patterns for right format
    gapped_patterns_sanitized = gapped_patterns.split()
    # print(f"{gapped_patterns_sanitized=}")

    # Create first and second patterns
    first_patterns = [pair.split('|')[0] for pair in gapped_patterns_sanitized]
    second_patterns = [pair.split('|')[1] for pair in gapped_patterns_sanitized]
    # print(f"\n\n{first_patterns=}")
    # print(f"\n\n{second_patterns=}")

    # Create prefix string
    # NOTE: The way this is done is that it takes the full first pattern and then adds the last letter from each subsequent pattern
    # You can also take just the first letter from each pattern and then put in the last full pattern
    # This is because suffix of first pattern = prefix of second pattern and so on
    # For example:
    # GACC, ACCG, CCGA, CGAG, GAGC ---> GACC G A G C
    # OR
    # GACC, ACCG, CCGA, CGAG, GAGC ---> G A C C GAGC
    prefix_string = first_patterns[0]
    # print(f"\n\n{prefix_string=}")
    for pattern in first_patterns[1:]:
        prefix_string += pattern[-1]
        # print(f"\n\n{pattern=}")
        # print(f"{pattern[-1]=}")
        # print(f"{prefix_string=}")

    # Create suffix string
    suffix_string = second_patterns[0]
    for pattern in second_patterns[1:]:
        suffix_string += pattern[-1]
        # print(f"\n\n{pattern=}")
        # print(f"{pattern[-1]=}")
        # print(f"{suffix_string=}")

    # Check for perfect overlap
    for i in range(k + d, len(prefix_string)):
        if prefix_string[i] != suffix_string[i - k - d]:
            return "there is no string spelled by the gapped patterns"

    # Return final answer
    return prefix_string + suffix_string[-(k + d):]


if __name__ == "__main__":

# # Sample dataset
#     k = 4
#     d = 2
#     gapped_patterns = 'GACC|GCGC ACCG|CGCC CCGA|GCCG CGAG|CCGG GAGC|CGGA'
#     # gapped_patterns_2 = 'GAGA|TTGA TCGT|GATG CGTG|ATGT TGGT|TGAG GTGA|TGTT GTGG|GTGA TGAG|GTTG GGTC|GAGA GTCG|AGAT'
#     answer = string_spelled_by_gapped_patterns(gapped_patterns, k, d)
#     print("\n", answer)

# # From file
    # Get dataset
    from pathlib import Path as partho
    current_dir = partho(__file__).parent
    filename = input("Please enter the filename: ")
    file_path = current_dir / filename

    with open(file_path, 'r') as file:
       # Read kmer length
       k, d = file.readline().split()
       gapped_patterns = file.readline().strip()

    chullu = string_spelled_by_gapped_patterns(gapped_patterns, int(k), int(d))
    print(chullu)
    # The above is done to close the loop and make it circular
    # as the overlap is k-1, so you leave out k-1 nucleotides at the end


    with open("Wk2_6_output.txt", "w") as output_file:
        output_file.write(chullu)

