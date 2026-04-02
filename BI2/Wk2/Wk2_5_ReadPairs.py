#########################################################################################
# Genome Sequencing  - Bioinformatic II Course from Coursera
#
# Week 2 - Composition from paired composition
# Given a string Text, a (k,d)-mer is a pair of k-mers in Text
# separated by distance d.
# We use the notation (Pattern1|Pattern2) to refer to a (k,d)-mer
# whose k-mers are Pattern1 and Pattern2.
# For example, (AAT|TGG) is a (3,4)-mer in TAATGCCATGGGATGTT.
# The (k,d)-mer composition of Text, denoted PairedCompositionk,d(Text),
# is the collection of all (k,d)- mers in Text (including repeated (k,d)-mers).
#########################################################################################
#########################################################################################

def read_pair_composition(genome, k, d):
    read_pairs = []
    for i in range(len(genome)-(2*k) -d+1):
        kmer1 = genome[i:i+k]
        kmer2 = genome[i+d+k:i+d+(2*k)]
        read_pairs.append((kmer1, kmer2))
    return read_pairs

answer = read_pair_composition('TAATGCCATGGGATGTT', 3, 2)
output_string = ' '.join([f"({k1}|{k2})" for k1, k2 in sorted(answer)])
print(output_string)