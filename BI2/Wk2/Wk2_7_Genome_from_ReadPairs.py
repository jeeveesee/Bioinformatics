#########################################################################################
# Genome Sequencing  - Bioinformatic II Course from Coursera
#
# Week 2 - Reconstruct string from read pairs where inputs are NOT in the correct order
#   Input: Integers k and d followed by a collection of paired k-mers PairedReads.
#   Output: A string Text with (k, d)-mer composition equal to PairedReads.
# To solve the String Reconstruction from Read-Pairs Problem,
# you will need to reconstruct a string from its path in the paired de Bruijn graph.
#########################################################################################

def genome_from_readpairs(k, d, read_pairs):
    """Recreates full genome from a list (or space-separated string) of read pairs.

    Implements reconstruction by building the paired de Bruijn graph
    and finding an Eulerian path. Handles inputs where read pairs are
    not in order.

    :param k: length of kmers
    :param d: gap between the paired kmers
    :param read_pairs: either a space-separated string of pairs "a|b ..."
                       or a list of strings ["a|b", ...]
    :return: reconstructed genome string
    """
    # Normalize input to list of "a|b" strings
    if isinstance(read_pairs, str):
        pairs = read_pairs.split()
    else:
        pairs = list(read_pairs)

    def build_paired_debruijn(pairs):
        """Return adjacency dict, indegree and outdegree maps for paired de Bruijn graph.

        Nodes are represented as strings 'p|q' where p and q are (k-1)-mers.
        """
        adj = {}
        indeg = {}
        outdeg = {}

        for pair in pairs:
            if '|' not in pair:
                raise ValueError(f"Invalid pair: {pair}")
            a, b = pair.split('|')
            src = a[:-1] + '|' + b[:-1]
            dst = a[1:] + '|' + b[1:]

            adj.setdefault(src, []).append(dst)

            outdeg[src] = outdeg.get(src, 0) + 1
            indeg[dst] = indeg.get(dst, 0) + 1

            # ensure keys exist for isolated nodes (explicit zero for clarity)
            indeg.setdefault(src, 0)
            outdeg.setdefault(dst, 0)

        return adj, indeg, outdeg

    def find_eulerian_path(adj, indeg, outdeg):
        """Find an Eulerian path in directed graph given adjacency and degree maps.

        Returns the list of nodes in path order.
        This is similar to what I have done. See the understanding text file to understand
        """
        # choose start node
        start = None
        for node in set(list(indeg.keys()) + list(outdeg.keys()) + list(adj.keys())):
            out = outdeg.get(node, 0)
            inn = indeg.get(node, 0)
            if out - inn == 1:
                start = node
                break

        if start is None:
            # pick any node with outgoing edges
            if adj:
                start = next(iter(adj))
            else:
                return []

        # make local copy of adjacency lists (we'll pop from them).
        # To get deterministic tie-breaking we sort neighbours lexicographically.
        # Because we use pop() (which removes from the end), we sort in
        # reverse order so that the smallest neighbor is popped last and
        # therefore consumed first during traversal. This yields a
        # lexicographically smallest Eulerian traversal when multiple edges
        # are available from a node.
        local_adj = {u: sorted(vs, reverse=True) for u, vs in adj.items()}
        # print("\n")
        # print(f"{local_adj=}")
        # print(f"{start=}")

        stack = [start]
        # print(f"{stack=}")
        path = []

        while stack:
            v = stack[-1]
            # print(f"{v=}")
            if local_adj.get(v):
                u = local_adj[v].pop()
                # print(f"{u=}")
                stack.append(u)
                # print(f"{stack=}")
            else:
                path.append(stack.pop())
                # print(f"Partho{path=}")
            # print(f"{local_adj=}")

        path.reverse()
        return path

    adj, indeg, outdeg = build_paired_debruijn(pairs)
    path = find_eulerian_path(adj, indeg, outdeg)
    # print(f"{path=}")

    if not path:
        return ""

    # Reconstruct prefix and suffix strings from path nodes
    prefix_nodes = [node.split('|')[0] for node in path]
    suffix_nodes = [node.split('|')[1] for node in path]
    # print(f"{prefix_nodes=}")
    # print(f"{suffix_nodes=}")

    prefix_string = prefix_nodes[0] + ''.join(p[-1] for p in prefix_nodes[1:])
    suffix_string = suffix_nodes[0] + ''.join(s[-1] for s in suffix_nodes[1:])
    # print(f"{prefix_string=}")
    # print(f"{suffix_string=}")

    # The final genome is prefix_string + last (k+d) chars of suffix_string
    if len(suffix_string) < k + d:
        raise ValueError("Suffix string too short to assemble with given k and d")

    # print(f"{suffix_string[-(k+d):]=}")
    genome = prefix_string + suffix_string[-(k + d):]
    return genome

if __name__ == "__main__":

# Sample dataset
    # k = 4
    # d = 2
    # read_pairs = 'GAGA|TTGA TCGT|GATG CGTG|ATGT TGGT|TGAG GTGA|TGTT GTGG|GTGA TGAG|GTTG GGTC|GAGA GTCG|AGAT'
    # gapped_patterns_2 = 'GAGA|TTGA TCGT|GATG CGTG|ATGT TGGT|TGAG GTGA|TGTT GTGG|GTGA TGAG|GTTG GGTC|GAGA GTCG|AGAT'
    k_exam = 3
    d_exam = 1
    read_pairs_exam = 'ACC|ATA ACT|ATT ATA|TGA ATT|TGA CAC|GAT CCG|TAC CGA|ACT CTG|AGC CTG|TTC GAA|CTT GAT|CTG GAT|CTG TAC|GAT TCT|AAG TGA|GCT TGA|TCT TTC|GAA'
    # Expected answer = CACCGATACTGATTCTGAAGCTT
    answer = genome_from_readpairs(k_exam, d_exam, read_pairs_exam)
    print("\n", answer)

# # # From file
#     # Get dataset
#     from pathlib import Path as partho
#     current_dir = partho(__file__).parent
#     filename = input("Please enter the filename: ")
#     file_path = current_dir / filename

#     with open(file_path, 'r') as file:
#        # Read kmer length
#        k, d = file.readline().split()
#        read_pairs = file.readline().strip()

#     chullu = genome_from_readpairs(int(k), int(d), read_pairs)
#     print(chullu)
#     # The above is done to close the loop and make it circular
#     # as the overlap is k-1, so you leave out k-1 nucleotides at the end


#     with open("Wk2_7_output.txt", "w") as output_file:
#         output_file.write(chullu)
