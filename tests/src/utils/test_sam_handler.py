from DAJIN2.utils.sam_handler import split_cigar
from DAJIN2.utils.sam_handler import calculate_alignment_length


def test_split_cigar():
    assert split_cigar("10M2D3M") == ["10M", "2D", "3M"]
    assert split_cigar("5M1I5M") == ["5M", "1I", "5M"]
    assert split_cigar("3M1D1I3M") == ["3M", "1D", "1I", "3M"]
    assert split_cigar("10M") == ["10M"]
    assert split_cigar("") == []
    assert split_cigar("3S10M2D3M") == ["3S", "10M", "2D", "3M"]
    assert split_cigar("5M1I5M3H") == ["5M", "1I", "5M", "3H"]
    assert split_cigar("3S3M1D1I3M3H") == ["3S", "3M", "1D", "1I", "3M", "3H"]


def test_calculate_alignment_length():
    assert calculate_alignment_length("10M2D3M") == 15
    assert calculate_alignment_length("5M1I5M") == 10
    assert calculate_alignment_length("3M1D1I3M") == 7
    assert calculate_alignment_length("3S10M2D3M") == 15
    assert calculate_alignment_length("5M1I5M3H") == 10
    assert calculate_alignment_length("3S3M1D1I3M3H") == 7
    assert calculate_alignment_length("10M") == 10
    assert calculate_alignment_length("") == 0
