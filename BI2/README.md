# Genome from Read-Pairs (Bioinformatics II)

Short implementation for the "String Reconstruction from Read-Pairs" problem (paired de Bruijn graph + Eulerian path).

What this repo contains

- `Wk2/Wk2_7_Genome_from_ReadPairs.py` — implementation of the assembler. Function `genome_from_readpairs(k, d, read_pairs)` accepts either a space-separated string of `a|b` pairs or a list of such strings.
- `tests/test_genome_from_readpairs.py` — pytest suite (5 tests) covering sample input, list input, malformed input, short-suffix error, and shuffled-input invariance.

Quick usage

Run the example script (prints assembled genome for the included sample):

```powershell
python "Wk2\Wk2_7_Genome_from_ReadPairs.py"
```

From the repository root you can also import the function in your own script:

```python
from Wk2.Wk2_7_Genome_from_ReadPairs import genome_from_readpairs

genome = genome_from_readpairs(k, d, read_pairs)
```

Run tests

Requires Python and pytest. From the project root:

```powershell
python -m pytest -q
```

Notes and next steps

- The adjacency ordering is currently non-deterministic; if you need strictly reproducible results across Python versions/environments, sort adjacency lists before traversal.
- I can add a small CLI to read `k`, `d`, and read-pairs from a file and write output to a file if you'd like.

Contact / author

- Files were edited in this workspace to implement the solution and tests; open an issue or request further changes if needed.
