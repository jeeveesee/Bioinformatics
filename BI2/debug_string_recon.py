from Wk1.Wk1_6_deBruijnGraph_kmer import debruijn_graph_from_kmers
from Wk2.Wk2_2_EulerianPath import eulerian_path

kmers = 'AAAT AATG ACCC ACGC ATAC ATCA ATGC CAAA CACC CATA CATC CCAG CCCA CGCT CTCA GCAT GCTC TACG TCAC TCAT TGCA'
kmers_sanitized = kmers.strip().split()

dB = debruijn_graph_from_kmers(kmers_sanitized)
print('debruijn graph (dict):')
for k, v in sorted(dB.items()):
    print(f"{k}: {v}")

# Sanitize
adj = [f"{key}: {' '.join(value)}" for key, value in dB.items()]
print('\nadj list (sanitized):')
for line in sorted(adj):
    print(line)

path = eulerian_path(adj)
print('\nEulerian path:')
print(path)

nodes = path.split('->')
print('\nnodes:')
print(nodes)

# Show reconstructed genome using the same path_to_genome approach
from Wk1.Wk1_2_PathToGenome import path_to_genome
print('\npath_to_genome result:')
print(path_to_genome([nodes]))
