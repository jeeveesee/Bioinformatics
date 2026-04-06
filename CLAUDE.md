# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

# Coding instructions

-  Utilize the comments at the top of the file to determine what needs to be done in the code
-  In some cases, there might be some pseudocode available within the comments section at the top of the file for you to use
-  There is a sample test within the name=="main" section to test your code output. Ensure it matches the output exactly

# Coding Style

-  Use PEP8 style for class names, function names
-  Do not use camel case anywhere
-  In each function definition, please provide docstrings with explanation, parameters, return definitions
-  Make sure the sample test input under the name=="main" section is parsed properly for input into functions
-  Make sure the output is formatted to match the expected output under the sample test
- Pure Python, no external libraries (only `pathlib`, `collections`, `typing` from stdlib)
- Functions only — no classes
- Pseudocode from the textbook is quoted in the module docstring at the top of each file
- `from pathlib import Path as partho` is the local convention for `Path`
- Debug `print()` statements are commented out, not deleted
- Signed integers formatted with `f"{n:+d}"` when writing back to string form

## Repository Overview

Solutions to the three-course Bioinformatics specialization on Coursera (Compeau & Pevzner). Each course has its own top-level folder:

- **BI1** — DNA pattern matching, motif finding (Wk_1–Wk_5)
- **BI2** — Genome assembly, de Bruijn graphs, Eulerian paths (Wk1–Wk4)
- **BI3** — Sequence alignment, dynamic programming, genome rearrangements (Wk1–Wk4)

## Running Code

Scripts are designed as standalone executables. Most have an embedded sample test in the `if __name__ == "__main__":` block and a commented-out file-reading section for dataset input.

```bash
# BI1 — run directly from the repo root
python BI1/Wk_1/Wk1_1_hidden_messages.py

# BI2 — run directly or as a module from the BI2 folder
python BI2/Wk2/Wk2_7_Genome_from_ReadPairs.py

# BI3 — MUST use module syntax from the BI3 folder (due to sibling imports)
cd BI3
python -m Wk4.Wk4_9_2BreakDistance   # no .py extension
```


## Architecture

### Module dependencies

Some weeks build on prior weeks via direct sibling imports. Each file adds `sys.path.insert(0, str(Path(__file__).parent))` to enable this, for example:

```
Wk4_9_2BreakDistance
  └── Wk4_8_2BreakOnGenome
        ├── Wk4_5_ColoredEdges → Wk4_3_Chromo2Cycle
        ├── Wk4_6_Graph2Genome → Wk4_4_Cycle2Chromo
        └── Wk4_7_2BreakOnGenomeGraph
```

This means such scripts (e.g., if one exists in BI3 folder) **must** be invoked with `python -m` from the `BI3/` directory, not run directly.

### File I/O convention

- Input files: `WkX_Y_input.txt` alongside the script; read via `current_dir = Path(__file__).parent`
- Output files: `WkX_Y_output.txt` written to the same directory
- Dataset input is always commented out in `__main__`; uncomment and call `input()` for the filename

### Genome representation

Throughout BI3 Wk4, chromosomes are represented in multiple forms that the modules translate between:

| Form | Example | Used by |
|------|---------|---------|
| Raw string | `"(+1 -2 +3)"` | file I/O, `__main__` |
| Cycle node list | `[1, 2, 4, 3, 5, 6]` | `Chromo2Cycle`, `ColoredEdges` |
| Colored edge list | `[(2,4), (3,6), (5,1)]` | `ColoredEdges`, `Graph2Genome` |
| Chromosome int list | `[1, -2, 3]` | `Cycle2Chromo`, `Graph2Genome` return |




