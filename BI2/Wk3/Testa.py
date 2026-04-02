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
    cyclic_spec = [0]
    n = len(peptide)
    
    # Generate all subpeptide masses (your original logic - which is correct!)
    for i in range(n):
        for j in range(i + 1, n + 1):
            cyclic_spec.append(prefix_mass[j] - prefix_mass[i])
            # Add wraparound subpeptide (complement) for internal segments
            if i > 0 and j < n:
                cyclic_spec.append(peptide_mass - (prefix_mass[j] - prefix_mass[i]))
    
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


###############################################################################################

if __name__ == "__main__":

# Sample test
    spectrum = [0, 113, 128, 186, 241, 299, 314, 427]
    # Expected answer = "186-128-113 186-113-128 128-186-113 128-113-186 113-186-128 113-128-186"
    answer = cyclopeptide_sequencing(spectrum)
    ## print(''.join(answer))
    print(answer)