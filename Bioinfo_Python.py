"""
Branch-and-Bound Cyclopeptide Sequencing Algorithm in Python
With performance benchmarking and comparison features
"""

import time
from collections import Counter
from typing import List, Tuple

# Amino acid mass table (integer masses)
AMINO_ACID_MASS = {
    'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103,
    'I': 113, 'L': 113, 'N': 114, 'D': 115, 'K': 128, 'Q': 128,
    'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186
}

def mass(peptide: Tuple[int, ...]) -> int:
    """Calculate the mass of a peptide (tuple of amino acid masses)"""
    return sum(peptide)

def parent_mass(spectrum: List[int]) -> int:
    """Get parent mass (largest value in spectrum)"""
    return max(spectrum)

def linear_spectrum(peptide: Tuple[int, ...]) -> List[int]:
    """Generate theoretical linear spectrum"""
    if len(peptide) == 0:
        return [0]
    
    prefix_mass = [0]
    for mass_val in peptide:
        prefix_mass.append(prefix_mass[-1] + mass_val)
    
    linear_spec = [0]
    for i in range(len(peptide)):
        for j in range(i + 1, len(peptide) + 1):
            linear_spec.append(prefix_mass[j] - prefix_mass[i])
    
    return sorted(linear_spec)

def cyclospectrum(peptide: Tuple[int, ...]) -> List[int]:
    """Generate theoretical cyclic spectrum"""
    if len(peptide) == 0:
        return [0]
    
    prefix_mass = [0]
    for mass_val in peptide:
        prefix_mass.append(prefix_mass[-1] + mass_val)
    
    peptide_mass = prefix_mass[-1]
    cyclic_spec = [0, peptide_mass]
    
    for i in range(len(peptide)):
        for j in range(i + 1, len(peptide) + 1):
            subpeptide_mass = prefix_mass[j] - prefix_mass[i]
            cyclic_spec.append(subpeptide_mass)
            if i > 0 and j < len(peptide):
                cyclic_spec.append(peptide_mass - subpeptide_mass)
    
    return sorted(cyclic_spec)

def is_consistent(peptide: Tuple[int, ...], spectrum: List[int]) -> bool:
    """Check if peptide is consistent with spectrum"""
    linear_spec = linear_spectrum(peptide)
    
    # Create frequency counters
    spectrum_freq = Counter(spectrum)
    linear_freq = Counter(linear_spec)
    
    # Check if all masses in linear spectrum are in experimental spectrum
    # with sufficient frequency
    for mass_val, count in linear_freq.items():
        if spectrum_freq[mass_val] < count:
            return False
    
    return True

def expand(peptides: List[Tuple[int, ...]], 
           masses: dict = AMINO_ACID_MASS) -> List[Tuple[int, ...]]:
    """Expand candidate peptides by adding each possible amino acid mass"""
    if len(peptides) == 0:
        return [tuple([m]) for m in set(masses.values())]
    
    expanded = []
    unique_masses = set(masses.values())
    
    for peptide in peptides:
        for mass_val in unique_masses:
            expanded.append(peptide + (mass_val,))
    
    return expanded

def peptide_to_string(peptide: Tuple[int, ...]) -> str:
    """Convert peptide (mass tuple) to string representation"""
    return '-'.join(map(str, peptide))

def cyclopeptide_sequencing(spectrum: List[int]) -> List[str]:
    """
    Main cyclopeptide sequencing algorithm
    Returns list of peptide strings that match the spectrum
    """
    candidate_peptides = [tuple()]  # Start with empty peptide
    final_peptides = []
    parent = parent_mass(spectrum)
    
    while candidate_peptides:
        # Expand all candidate peptides
        candidate_peptides = expand(candidate_peptides)
        
        # Filter peptides
        new_candidates = []
        for peptide in candidate_peptides:
            peptide_mass = mass(peptide)
            
            if peptide_mass == parent:
                # Check if cyclospectrum matches
                if cyclospectrum(peptide) == spectrum:
                    peptide_str = peptide_to_string(peptide)
                    if peptide_str not in final_peptides:
                        final_peptides.append(peptide_str)
                # Don't add to new_candidates (remove from search)
            elif is_consistent(peptide, spectrum):
                # Keep peptide if it's consistent
                new_candidates.append(peptide)
            # else: remove from search (not consistent)
        
        candidate_peptides = new_candidates
    
    return final_peptides

def benchmark_algorithm(spectrum: List[int], runs: int = 5) -> dict:
    """Benchmark the algorithm and return statistics"""
    times = []
    
    for _ in range(runs):
        start = time.perf_counter()
        result = cyclopeptide_sequencing(spectrum)
        end = time.perf_counter()
        times.append(end - start)
    
    return {
        'result': result,
        'avg_time': sum(times) / len(times),
        'min_time': min(times),
        'max_time': max(times),
        'times': times
    }

# Example usage and benchmarking
if __name__ == "__main__":
    print("=" * 60)
    print("Cyclopeptide Sequencing - Python Implementation")
    print("=" * 60)
    
    # Example 1: NQEL spectrum
    print("\nExample 1: Testing with spectrum corresponding to NQEL")
    example_spectrum = [0, 113, 114, 128, 129, 242, 242, 257, 370, 371, 484]
    
    print(f"Spectrum: {example_spectrum}")
    benchmark1 = benchmark_algorithm(example_spectrum, runs=5)
    
    print(f"\nFound peptides: {benchmark1['result']}")
    print(f"Average time: {benchmark1['avg_time']*1000:.4f} ms")
    print(f"Min time: {benchmark1['min_time']*1000:.4f} ms")
    print(f"Max time: {benchmark1['max_time']*1000:.4f} ms")
    
    # Example 2: Simple spectrum
    print("\n" + "=" * 60)
    print("\nExample 2: Testing with simple spectrum")
    test_spectrum = [0, 71, 87, 158, 158, 229, 229, 316]
    
    print(f"Spectrum: {test_spectrum}")
    benchmark2 = benchmark_algorithm(test_spectrum, runs=5)
    
    print(f"\nFound peptides: {benchmark2['result']}")
    print(f"Average time: {benchmark2['avg_time']*1000:.4f} ms")
    print(f"Min time: {benchmark2['min_time']*1000:.4f} ms")
    print(f"Max time: {benchmark2['max_time']*1000:.4f} ms")
    
    # Example 3: Larger spectrum for performance testing
    print("\n" + "=" * 60)
    print("\nExample 3: Larger spectrum (performance test)")
    large_spectrum = [0, 71, 99, 101, 103, 128, 129, 170, 200, 202, 
                      228, 231, 257, 299, 303, 328, 330, 332, 333, 
                      403, 431, 434, 461, 534]
    
    print(f"Spectrum length: {len(large_spectrum)}")
    print(f"Parent mass: {max(large_spectrum)}")
    
    benchmark3 = benchmark_algorithm(large_spectrum, runs=3)
    
    print(f"\nFound {len(benchmark3['result'])} peptide(s)")
    if len(benchmark3['result']) <= 5:
        print(f"Peptides: {benchmark3['result']}")
    print(f"Average time: {benchmark3['avg_time']*1000:.4f} ms")
    print(f"Min time: {benchmark3['min_time']*1000:.4f} ms")
    print(f"Max time: {benchmark3['max_time']*1000:.4f} ms")
    
    print("\n" + "=" * 60)
    print("\nPerformance Comparison Notes:")
    print("- Python uses tuples (immutable) for efficient hashing")
    print("- Counter class provides O(1) frequency lookups")
    print("- List comprehensions offer fast iteration")
    print("- Run this alongside the R implementation to compare!")
    print("=" * 60)