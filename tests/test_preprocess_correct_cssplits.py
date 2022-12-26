"""
- repetitive dels
- Shared mutation percentages in 5-mer
    - sample:  [{"=A,=A,=A,=A,=A": [0.98,0.98,0.98,0.98,0.98,0.98]}, ]
    - control: [{"=A,=A,=A,=A,=A": [0.98,0.98,0.98,0.98,0.98,0.98]}]
"""

from misc.scripts import scratch_correct_cssplits


def test_extract_indexes_with_both_ends_not_N___noN():
    cssplits = "=A,=C,=G,=T"
    test = scratch_correct_cssplits.extract_indexes_with_both_ends_not_N(cssplits)
    answer = (0, 3)
    assert test == answer


def test_extract_indexes_with_both_ends_not_N___leftN():
    cssplits = "N,N,N,N,=A,=C,=G,=T"
    test = scratch_correct_cssplits.extract_indexes_with_both_ends_not_N(cssplits)
    answer = (4, 7)
    assert test == answer


def test_extract_indexes_with_both_ends_not_N___rightN():
    cssplits = "=A,=C,=G,=T,N,N,N,N"
    test = scratch_correct_cssplits.extract_indexes_with_both_ends_not_N(cssplits)
    answer = (0, 3)
    assert test == answer


def test_extract_indexes_with_both_ends_not_N___bothN():
    cssplits = "N,N,=A,=C,=G,=T,=A,N,N"
    test = scratch_correct_cssplits.extract_indexes_with_both_ends_not_N(cssplits)
    answer = (2, 6)
    assert test == answer
