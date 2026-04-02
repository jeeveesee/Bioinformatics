#########################################################################################
# Genome Sequencing  - Bioinformatic III Course from Coursera
#
# Week 4 - Random Breakage Model of Chromosome Evolution - Number of Breakpoints
# Code Challenge: Find the number of breakpoints in a permutation.
# Input: A permutation.
# Output: The number of breakpoints in this permutation.
#########################################################################################


def number_of_breakpoints(permutation):
    """
    Counts the number of breakpoints in a signed permutation.
    A breakpoint is any adjacent pair in the extended permutation (0, p1..pn, n+1)
    where the next element is not exactly 1 more than the current element.

    Parameters:
    permutation -> list of strings with +/- signs e.g. ["+3", "+4", "+5", "-12"]

    Returns:
    integer count of breakpoints
    """
    # Convert to integers and build the extended permutation with 0 and n+1 bookends
    perm = [int(x) for x in permutation]
    n = len(perm)
    extended = [0] + perm + [n + 1]

    # Count adjacent pairs where the next element is not exactly current + 1
    breakpoints = 0
    for i in range(len(extended) - 1):
        if extended[i + 1] - extended[i] != 1:  # not a proper adjacency
            breakpoints += 1

    return breakpoints


###########################################################################

if __name__ == "__main__":
    # Sample test
    # permutation = [
    #     "+3",
    #     "+4",
    #     "+5",
    #     "-12",
    #     "-8",
    #     "-7",
    #     "-6",
    #     "+1",
    #     "+2",
    #     "+10",
    #     "+9",
    #     "-11",
    #     "+13",
    #     "+14",
    # ]
    # # Expected answer = 8
    # answer = number_of_breakpoints(permutation)
    # print(answer)

    # # From file

    # # Get dataset
    # from pathlib import Path as partho

    # current_dir = partho(__file__).parent
    # filename = input("Please enter the filename: ")
    # file_path = current_dir / filename

    # with open(file_path, "r") as file:
    #     permutation = file.read().strip().split()

    # answer = number_of_breakpoints(permutation)
    # print(answer)

    # with open("Wk4_2_output.txt", "w") as output_file:
    #     output_file.write(str(answer))

    # Exam
    permutation = ["+6", "-12", "-9", "+17", "+18", "-4", "+5",
                   "-3", "+11", "+19", "+20", "+10", "+8", "+15", "-14",
                    "-13", "+2", "+7", "-16", "-1"]
    # Expected answer = 8
    answer = number_of_breakpoints(permutation)
    print(answer)
