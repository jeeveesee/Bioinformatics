import pytest

from Wk2.Wk2_7_Genome_from_ReadPairs import genome_from_readpairs


def test_sample():
    k = 4
    d = 2
    read_pairs = (
        'GAGA|TTGA TCGT|GATG CGTG|ATGT TGGT|TGAG '
        'GTGA|TGTT GTGG|GTGA TGAG|GTTG GGTC|GAGA GTCG|AGAT'
    )
    expected = "GTGGTCGTGAGATGTTGA"
    assert genome_from_readpairs(k, d, read_pairs) == expected


def test_list_input_same_result():
    k = 4
    d = 2
    read_pairs = (
        'GAGA|TTGA TCGT|GATG CGTG|ATGT TGGT|TGAG '
        'GTGA|TGTT GTGG|GTGA TGAG|GTTG GGTC|GAGA GTCG|AGAT'
    )
    pairs_list = read_pairs.split()
    expected = "GTGGTCGTGAGATGTTGA"
    assert genome_from_readpairs(k, d, pairs_list) == expected


def test_invalid_pair_raises():
    k = 3
    d = 1
    bad = 'AAA BBB|CCC'
    with pytest.raises(ValueError):
        genome_from_readpairs(k, d, bad)


def test_suffix_too_short_raises():
    # Single pair will produce suffix shorter than k+d and should raise
    k = 3
    d = 1
    single = 'AAA|AAA'
    with pytest.raises(ValueError):
        genome_from_readpairs(k, d, single)


def test_shuffled_input_same_result():
    k = 4
    d = 2
    read_pairs = (
        'GAGA|TTGA TCGT|GATG CGTG|ATGT TGGT|TGAG '
        'GTGA|TGTT GTGG|GTGA TGAG|GTTG GGTC|GAGA GTCG|AGAT'
    )
    pairs_list = read_pairs.split()
    shuffled = list(reversed(pairs_list))
    expected = "GTGGTCGTGAGATGTTGA"
    assert genome_from_readpairs(k, d, shuffled) == expected
