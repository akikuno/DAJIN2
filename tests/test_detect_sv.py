from importlib import reload
import pytest
from src.DAJIN2 import classification

reload(classification)


# def test_match():
#     test = "=A,=C,=G,=T,=A,=A,=A,=A,=T,=G,=C,=A"
#     test = detect_sv.detect_sv(test, kmer_size=4, window_size=2, threshold=2)
#     answer = False
#     assert test == answer


# def test_snv():
#     test = "=A,=C,=G,=T,*AT,=A,=A,=A,=T,=G,=C,=A"
#     test = detect_sv.detect_sv(test, kmer_size=4, window_size=2, threshold=3)
#     answer = False
#     assert test == answer


# def test_small_indel():
#     test = "=A,=C,=G,=T,-A,+N|=A,=A,=T,=G,=C,=A"
#     test = detect_sv.detect_sv(test, kmer_size=4, window_size=2, threshold=3)
#     answer = False
#     assert test == answer


# def test_large_deletion():
#     test = "=A,=C,=G,=T,-A,-A,=A,=A,=T,=G,=C,=A"
#     test = detect_sv.detect_sv(test, kmer_size=4, window_size=2, threshold=2)
#     answer = True
#     assert test == answer


# def test_valueerror():
#     with pytest.raises(ValueError) as excinfo:
#         test = "=A,=C,=G,=T,-A,-A,=A,=A,=T,=G,=C,=A"
#         test = detect_sv.detect_sv(test, kmer_size=1, window_size=100, threshold=2)
#     assert "kmer_size must be larger than window_size" in str(excinfo.value)


def test_detect_sv_N():
    test = "N,N,=A,=A"
    test = classification.detect_sv(test, threshold=2)
    answer = True
    assert test == answer


def test_detect_sv_insertion_true():
    test = "=A,=A,+A|+A|=A"
    test = classification.detect_sv(test, threshold=2)
    answer = True
    assert test == answer


def test_detect_sv_insertion_false():
    test = "=A,=A,+A|+A|=A"
    test = classification.detect_sv(test, threshold=3)
    answer = False
    assert test == answer


def test_detect_sv_deletion():
    test = "=A,=A,-A,-A"
    test = classification.detect_sv(test, threshold=2)
    answer = True
    assert test == answer


def test_detect_sv_substitution():
    test = "=A,=A,*AG,*AG"
    test = classification.detect_sv(test, threshold=2)
    answer = True
    assert test == answer

