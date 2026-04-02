def spectrum_score(theoretical_spectrum: list[int], experimental_spectrum: list[int]) -> int:
    """
    Compute the  score of a peptide represented as list of amino-acid masses (can be linear or cyclical)

    Parameters:
        theoretical_spectrum -> list of integer representing theoretical spectrum masses (e.g. [113,128,186])
        experimental_spectrum -> list of integers representing experimental spectrum

    Returns:
        integer score (counts multiplicities)
    """
    
    theo_count = Counter(theoretical_spectrum)
    exper_count = Counter(experimental_spectrum)
    score = 0
    for mass in theo_count:
        if mass in exper_count:
            score += min(theo_count[mass], exper_count[mass])
    return score

print(spectrum_score([131, 71, 131, 71], [0, 71, 71, 71, 131, 131, 131, 156, 198, 199, 199, 202, 202, 202, 333, 333, 333, 404, 404]))