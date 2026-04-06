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
from BI3.Wk5.Wk5_3_ColoredEdges import colored_edges
from BI3.Wk5.Wk5_6_2BreakOnGenome import two_break_on_genome


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
    # Sample test
    genome_P_raw = "(+1 -2 -3 +4)"
    genome_Q_raw = "(+1 +2 -4 -3)"

    # Expected answer =
    # (+1 -2 -3 +4)
    # (+1 -2 -3)(+4)
    # (+1 -2 -4 -3)
    # (-3 +1 +2 -4)
    answer = two_break_sorting(genome_P_raw, genome_Q_raw)
    print(answer)

    # # From file
    # # Get dataset
    # from pathlib import Path as partho

    # current_dir = partho(__file__).parent
    # filename = input("Please enter the filename: ")
    # file_path = current_dir / filename

    # with open(file_path, "r") as file:
    #     genome_P_raw = file.readline().strip()
    #     genome_Q_raw = file.readline().strip()
    #     # i1, i2, i3, i4 = map(int, file.readline().strip().split(","))

    # answer = two_break_distance(genome_P_raw, genome_Q_raw)
    # answer = [tuple(f"{n:+d}" for n in chrom) for chrom in answer]
    # answer = "".join("(" + " ".join(chrom) + ")" for chrom in answer)
    # print(answer)

    # # answer = "(" + " ".join(map(str, answer)) + ")"
    # # answer = "(" + " ".join(f"{n:+d}" for n in answer) + ")"
    # # answer = ", ".join(str(e) for e in answer)
    # # answer = "".join("(" + " ".join(f"{n:+d}" for n in chrom) + ")" for chrom in answer)

    # with open(current_dir / "Wk4_10_output.txt", "w") as output_file:
    #     output_file.write(str(answer))


    # Exam
    # permutation = ["+6", "-12", "-9", "+17", "+18", "-4", "+5",
    #                "-3", "+11", "+19", "+20", "+10", "+8", "+15", "-14",
    #                 "-13", "+2", "+7", "-16", "-1"]
    # # Expected answer = 8
    # answer = number_of_breakpoints(permutation)
    # print(answer)



"""
Notes:

Step 1 — Get colored edges for P and Q

  colored_edges() (from Wk4_5) converts each genome into a list of node pairs. Each gene i maps to two nodes:   2i-1 (tail) and 2i (head). So:

  ┌──────┬───────────┬───────────┐
  │ Gene │ Tail node │ Head node │
  ├──────┼───────────┼───────────┤
  │ +1   │ 1         │ 2         │
  ├──────┼───────────┼───────────┤
  │ +2   │ 3         │ 4         │
  ├──────┼───────────┼───────────┤
  │ +3   │ 5         │ 6         │
  ├──────┼───────────┼───────────┤
  │ +4   │ 7         │ 8         │
  ├──────┼───────────┼───────────┤
  │ +5   │ 9         │ 10        │
  ├──────┼───────────┼───────────┤
  │ +6   │ 11        │ 12        │
  └──────┴───────────┴───────────┘

  Colored edges connect the head of gene j to the tail of gene j+1 (wrapping around):

  edges_P = [(2,3), (4,5), (6,7), (8,9), (10,11), (12,1)]
  edges_Q = [(2,6), (5,12), (11,10), (9,1), (4,8), (7,3)]

  ---
  Step 2 — Build adjacency maps

  We turn those edge lists into dictionaries so we can say "from node X, where do I go?":

  adj_P:  2↔3,  4↔5,  6↔7,  8↔9,  10↔11,  12↔1
  adj_Q:  2↔6,  5↔12, 11↔10, 9↔1,  4↔8,   7↔3

  ---
  Step 3 — Count cycles by alternating P→Q

  In the breakpoint graph, edges from P and Q alternate to form cycles. To find a cycle, we:
  1. Start at any unvisited node
  2. Follow a P-edge
  3. Follow a Q-edge
  4. Repeat until we return to the start

  Cycle 1 — start at node 2:
  2 --[P]--> 3 --[Q]--> 7 --[P]--> 6 --[Q]--> 2  ✓ back to start
  Nodes visited: {2, 3, 6, 7}

  Cycle 2 — next unvisited is 4:
  4 --[P]--> 5 --[Q]--> 12 --[P]--> 1 --[Q]--> 9 --[P]--> 8 --[Q]--> 4  ✓
  Nodes visited: {4, 5, 8, 9, 1, 12}

  Cycle 3 — next unvisited is 10:
  10 --[P]--> 11 --[Q]--> 10  ✓ back to start immediately
  Nodes visited: {10, 11}

  Total cycles = 3

  ---
  Step 4 — Final answer

  d(P, Q) = blocks − cycles = 6 − 3 = 3  ✓

  ---
  Why does this formula work?

  Intuitively: if P and Q were identical, every gene would form its own 2-node cycle → cycles = blocks →     
  distance = 0. Every 2-break operation can increase the cycle count by at most 1, so the minimum number of  
  2-breaks needed to transform P into Q is exactly how many cycles are "missing" from the maximum possible. 

"""