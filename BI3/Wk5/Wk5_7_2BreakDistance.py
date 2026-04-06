#########################################################################################
# Genome Sequencing  - Bioinformatic III Course from Coursera
#
# Week 4 - 2-break distance
# Code Challenge: Solve the 2-Break Distance Problem.
# Input: Genomes P and Q.
# Output: The 2-break distance d(P, Q).
"""
2-BreakDistance(P, Q)
     return Blocks(P, Q) - Cycles(BreakpointGraph(P, Q))
"""
# NOTE: YOU HAVE TO RUN THIS FROM THE BI3 FOLDER WITH THE COMMAND:
# python -m Wk5.Wk5_7_2BreakDistance (NO .py)
#########################################################################################

# Main code

import sys
from pathlib import Path as partho

sys.path.insert(0, str(partho(__file__).parent))
from Wk5.Wk5_3_ColoredEdges import colored_edges


def two_break_distance(genome_p_raw, genome_q_raw):
    """
    Computes the 2-break distance between two genomes P and Q.

    Builds the breakpoint graph from the colored edges of both genomes,
    counts the number of alternating cycles, then returns:
        d(P, Q) = blocks(P, Q) - cycles(BreakpointGraph(P, Q))

    Parameters:
        genome_p_raw (str): Genome P as a string, e.g. "(+1 +2 +3 +4 +5 +6)".
        genome_q_raw (str): Genome Q as a string, e.g. "(+1 -3 -6 -5)(+2 -4)".

    Returns:
        int: The 2-break distance d(P, Q).
    """
    edges_p = colored_edges(genome_p_raw)
    edges_q = colored_edges(genome_q_raw)
    # print(f"{edges_p=}")
    # print(f"{edges_q=}")

    # Number of synteny blocks equals the number of colored edges in either genome
    # d(P, Q) = n(=blocks) - cycles
    n = len(edges_p)

    # Build adjacency maps for the breakpoint graph
    adj_p = {}
    for u, v in edges_p:
        adj_p[u] = v
        adj_p[v] = u

    adj_q = {}
    for u, v in edges_q:
        adj_q[u] = v
        adj_q[v] = u

    # print(f"{adj_p=}")
    # print(f"{adj_q=}")
    # Count cycles by alternating along P-edges then Q-edges
    visited = set()
    cycles = 0

    for start in adj_p:
        if start in visited:
            continue
        cycles += 1
        current = start
        while True:
            visited.add(current)
            next_node = adj_p[current]   # follow P-colored edge
            visited.add(next_node)
            current = adj_q[next_node]   # follow Q-colored edge
            if current == start:
                break
            # print(visited)

    return n - cycles


###########################################################################

if __name__ == "__main__":
    # # Sample test
    # genome_P_raw = "(+1 +2 +3 +4 +5 +6)"
    # genome_Q_raw = "(+1 -3 -6 -5)(+2 -4)"

    # # Expected answer = 3
    # answer = two_break_distance(genome_P_raw, genome_Q_raw)
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

    answer = two_break_distance(genome_P_raw, genome_Q_raw)
    # answer = [tuple(f"{n:+d}" for n in chrom) for chrom in answer]
    # answer = "".join("(" + " ".join(chrom) + ")" for chrom in answer)
    # print(answer)

    # # answer = "(" + " ".join(map(str, answer)) + ")"
    # # answer = "(" + " ".join(f"{n:+d}" for n in answer) + ")"
    # # answer = ", ".join(str(e) for e in answer)
    # # answer = "".join("(" + " ".join(f"{n:+d}" for n in chrom) + ")" for chrom in answer)

    with open(current_dir / "Wk5_7_output.txt", "w") as output_file:
        output_file.write(str(answer))


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