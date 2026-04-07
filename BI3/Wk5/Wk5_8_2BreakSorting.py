#########################################################################################
# Genome Sequencing  - Bioinformatic III Course from Coursera
#
# Week 4 - 2-break sorting
# The Breakpoint Theorem guarantees that there must always be a 2-break
# increasing the number of red-blue cycles in the breakpoint graph by 1.
# 2-Break Sorting Problem: Find a shortest transformation of one genome
# into another by 2-breaks.
# Code Challenge: Solve the 2-Break Sorting Problem.
# Input: Two genomes with circular chromosomes on the same set of synteny blocks.
# Output: The sequence of genomes resulting from applying a shortest sequence of
# 2-breaks transforming one genome into the other.
"""
ShortestRearrangementScenario(P, Q)
     output P
     RedEdges ← ColoredEdges(P)
     BlueEdges ← ColoredEdges(Q)
     BreakpointGraph ← the graph formed by RedEdges and BlueEdges
     while BreakpointGraph has a non-trivial cycle Cycle
          (i1 , i2 , i3 , i4 ) ← path starting at arbitrary blue edge in nontrivial red-blue cycle
          RedEdges ← RedEdges with edges (i1, i2) and (i3, i4) removed
          RedEdges ← RedEdges with edges (i1, i4) and (i2, i3) added
          BreakpointGraph ← the graph formed by RedEdges and BlueEdges
          P ← 2-BreakOnGenome(P, i1 , i2 , i4 , i3 )
          output P
"""
# NOTE: YOU MAY HAVE TO RUN THIS FROM THE BI3 FOLDER WITH THE COMMAND:
# python -m Wk5.Wk5_8_2BreakSorting (NO .py)
#########################################################################################

# Main code

import sys
from pathlib import Path as partho

sys.path.insert(0, str(partho(__file__).parent))
from Wk5.Wk5_3_ColoredEdges import colored_edges
from Wk5.Wk5_6_2BreakOnGenome import two_break_on_genome


def format_genome(chromosomes):
    """
    Converts a list of chromosomes to a genome string.

    Parameters:
        chromosomes (list of list of int): Each inner list is a chromosome of signed
                                           integers, e.g. [[1, -2, -3], [4]].

    Returns:
        str: Genome string, e.g. "(+1 -2 -3)(+4)".
    """
    return "".join("(" + " ".join(f"{n:+d}" for n in chrom) + ")" for chrom in chromosomes)


def find_non_trivial_cycle(adj_red, adj_blue):
    """
    Finds the node sequence of any non-trivial cycle in the breakpoint graph.

    Traverses cycles by alternating along red (P) then blue (Q) edges. A trivial
    cycle has exactly 2 nodes (the same pair of nodes is connected by both a red
    and a blue edge). A non-trivial cycle has 4 or more nodes.

    Parameters:
        adj_red (dict): Adjacency map for red (P) colored edges; maps each node to its
                        red neighbor.
        adj_blue (dict): Adjacency map for blue (Q) colored edges; maps each node to
                         its blue neighbor.

    Returns:
        list of int: Ordered node sequence of a non-trivial cycle, or None if all
                     cycles are trivial.
    """
    visited = set()
    for start in sorted(adj_red):
        if start in visited:
            continue
        cycle = []
        current = start
        while True:
            cycle.append(current)
            visited.add(current)
            next_node = adj_red[current]    # follow red edge
            visited.add(next_node)
            cycle.append(next_node)
            current = adj_blue[next_node]   # follow blue edge
            if current == start:
                break
        if len(cycle) > 2:
            return cycle
    return None


def two_break_sorting(genome_p_raw, genome_q_raw):
    """
    Finds a shortest sequence of 2-breaks transforming genome P into genome Q.

    Repeatedly locates a non-trivial red-blue cycle in the breakpoint graph and
    applies a 2-break to split it, increasing the cycle count by 1 each step,
    until P equals Q (all cycles are trivial).

    Parameters:
        genome_p_raw (str): Genome P as a string, e.g. "(+1 -2 -3 +4)".
        genome_q_raw (str): Genome Q as a string, e.g. "(+1 +2 -4 -3)".

    Returns:
        str: Newline-separated sequence of all intermediate genomes from P to Q
             (inclusive of both endpoints).
    """
    edges_blue = colored_edges(genome_q_raw)
    adj_blue = {}
    for u, v in edges_blue:
        adj_blue[u] = v
        adj_blue[v] = u

    result = [genome_p_raw]
    current_p_raw = genome_p_raw

    while True:
        edges_red = colored_edges(current_p_raw)
        adj_red = {}
        for u, v in edges_red:
            adj_red[u] = v
            adj_red[v] = u

        cycle = find_non_trivial_cycle(adj_red, adj_blue)
        if cycle is None:
            break

        # Take 4 consecutive cycle nodes: (i1,i2) = red edge, (i2,i3) = blue edge,
        # (i3,i4) = red edge. The 2-break removes red edges (i1,i2) and (i3,i4)
        # and adds (i1,i4) and (i2,i3), splitting the non-trivial cycle into two.
        i1, i2, i3, i4 = cycle[0], cycle[1], cycle[2], cycle[3]

        # Per pseudocode: P ← 2-BreakOnGenome(P, i1, i2, i4, i3)
        new_chromosomes = two_break_on_genome(current_p_raw, i1, i2, i4, i3)
        current_p_raw = format_genome(new_chromosomes)
        result.append(current_p_raw)

    return "\n".join(result)


###########################################################################

if __name__ == "__main__":
    # # Sample test
    # genome_P_raw = "(+1 -2 -3 +4)"
    # genome_Q_raw = "(+1 +2 -4 -3)"

    # # Expected answer =
    # # (+1 -2 -3 +4)
    # # (+1 -2 -3)(+4)
    # # (+1 -2 -4 -3)
    # # (-3 +1 +2 -4)
    # answer = two_break_sorting(genome_P_raw, genome_Q_raw)
    # print(answer)

    # From file
    # Get dataset
    from pathlib import Path as partho

    current_dir = partho(__file__).parent
    filename = input("Please enter the filename: ")
    file_path = current_dir / filename

    with open(file_path, "r") as file:
        genome_P_raw = file.readline().strip()
        genome_Q_raw = file.readline().strip()
        # i1, i2, i3, i4 = map(int, file.readline().strip().split(","))

    answer = two_break_sorting(genome_P_raw, genome_Q_raw)
    # answer = [tuple(f"{n:+d}" for n in chrom) for chrom in answer]
    # answer = "".join("(" + " ".join(chrom) + ")" for chrom in answer)
    # print(answer)

    # # answer = "(" + " ".join(map(str, answer)) + ")"
    # # answer = "(" + " ".join(f"{n:+d}" for n in answer) + ")"
    # # answer = ", ".join(str(e) for e in answer)
    # # answer = "".join("(" + " ".join(f"{n:+d}" for n in chrom) + ")" for chrom in answer)

    with open(current_dir / "Wk5_8_output.txt", "w") as output_file:
        output_file.write(str(answer))


    # Exam
    # permutation = ["+6", "-12", "-9", "+17", "+18", "-4", "+5",
    #                "-3", "+11", "+19", "+20", "+10", "+8", "+15", "-14",
    #                 "-13", "+2", "+7", "-16", "-1"]
    # # Expected answer = 8
    # answer = number_of_breakpoints(permutation)
    # print(answer)



"""
Notes on working:

Sure! Using the sample inputs:
  P = (+1 -2 -3 +4)
  Q = (+1 +2 -4 -3)

  ---
  Setup — Build the Blue Edges (Q, fixed forever)

  colored_edges(Q) gives us adj_blue. These never change — Q is our target.

  adj_blue:  2↔3,  4↔8,  7↔6,  5↔1

  ---
  Iteration 1

  Red edges from current P = (+1 -2 -3 +4):
  adj_red:  2↔4,  3↔6,  5↔7,  8↔1

  Find non-trivial cycle (start at smallest unvisited node = 1, alternate red→blue):
  1 -[R]→ 8 -[B]→ 4 -[R]→ 2 -[B]→ 3 -[R]→ 6 -[B]→ 7 -[R]→ 5 -[B]→ 1  ✓
  cycle = [1, 8, 4, 2, 3, 6, 7, 5] — length 8, non-trivial ✓

  Extract i1,i2,i3,i4 from first 4 nodes:
  i1=1, i2=8  ← red edge
        i2=8, i3=4  ← blue edge (8→4)
              i3=4, i4=2  ← red edge

  Apply 2-break two_break_on_genome(P, i1=1, i2=8, i4=2, i3=4):
  - Remove red edges (1,8) and (4,2)=(2,4)
  - Add red edges (1,2) and (8,4)=(4,8)

  New P → (+1)(+2 -4 +3)

  ---
  Iteration 2

  Red edges from (+1)(+2 -4 +3):
  adj_red:  2↔1,  4↔8,  7↔5,  6↔3

  Find non-trivial cycle (start at 1, alternate red→blue):
  1 -[R]→ 2 -[B]→ 3 -[R]→ 6 -[B]→ 7 -[R]→ 5 -[B]→ 1  ✓
  cycle = [1, 2, 3, 6, 7, 5] — length 6, non-trivial ✓

  Note: node 4 forms a trivial cycle 4-[R]→8-[B]→4 — red and blue agree on (4,8), so already "solved".

  Extract i1,i2,i3,i4:
  i1=1, i2=2  ← red edge
        i2=2, i3=3  ← blue edge (2→3)
              i3=3, i4=6  ← red edge

  Apply 2-break two_break_on_genome(P, 1, 2, 6, 3):
  - Remove red edges (1,2) and (3,6)
  - Add red edges (1,6) and (2,3)

  New P → (+1 +2 -4 +3)

  ---
  Iteration 3

  Red edges from (+1 +2 -4 +3):
  adj_red:  2↔3,  4↔8,  7↔5,  6↔1

  Find non-trivial cycle (start at 1):
  1 -[R]→ 6 -[B]→ 7 -[R]→ 5 -[B]→ 1  ✓
  cycle = [1, 6, 7, 5] — length 4, non-trivial ✓

  Note: nodes 2,3 and 4,8 are both trivial now.

  Extract i1,i2,i3,i4:
  i1=1, i2=6  ← red edge
        i2=6, i3=7  ← blue edge (6→7)
              i3=7, i4=5  ← red edge

  Apply 2-break two_break_on_genome(P, 1, 6, 5, 7):
  - Remove red edges (1,6) and (7,5)=(5,7)
  - Add red edges (1,5) and (6,7)

  New P → (+1 +2 -4 -3) = Q ✓

  ---
  Iteration 4 — Termination

  Red edges now match blue edges exactly:
  adj_red:  2↔3,  4↔8,  7↔6,  5↔1
  adj_blue: 2↔3,  4↔8,  7↔6,  5↔1

  Every cycle is trivial (length 2). find_non_trivial_cycle returns None → loop exits.

  ---
  Final Output

  (+1 -2 -3 +4)     ← initial P
  (+1)(+2 -4 +3)    ← after 2-break 1
  (+1 +2 -4 +3)     ← after 2-break 2
  (+1 +2 -4 -3)     ← after 2-break 3 = Q

  3 steps = d(P, Q) = 3 ✓ — guaranteed shortest because each 2-break increases the cycle count by exactly 1, and we need blocks − cycles = 4 − 1 = 3 additional
  cycles to reach the trivial state. 

"""

