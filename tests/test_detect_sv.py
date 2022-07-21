from importlib import reload
import pytest
from src.DAJIN2.classification import detect_sv

reload(detect_sv)


def test_match():
    test = "=A,=C,=G,=T,=A,=A,=A,=A,=T,=G,=C,=A"
    test = detect_sv.is_sv(test, kmer_size=4, window_size=2, threshold=2)
    answer = False
    assert test == answer


def test_snv():
    test = "=A,=C,=G,=T,*AT,=A,=A,=A,=T,=G,=C,=A"
    test = detect_sv.is_sv(test, kmer_size=4, window_size=2, threshold=3)
    answer = False
    assert test == answer


def test_small_indel():
    test = "=A,=C,=G,=T,-A,+N|=A,=A,=T,=G,=C,=A"
    test = detect_sv.is_sv(test, kmer_size=4, window_size=2, threshold=3)
    answer = False
    assert test == answer


def test_large_deletion():
    test = "=A,=C,=G,=T,-A,-A,=A,=A,=T,=G,=C,=A"
    test = detect_sv.is_sv(test, kmer_size=4, window_size=2, threshold=2)
    answer = True
    assert test == answer


def test_valueerror():
    with pytest.raises(ValueError) as excinfo:
        test = "=A,=C,=G,=T,-A,-A,=A,=A,=T,=G,=C,=A"
        test = detect_sv.is_sv(test, kmer_size=1, window_size=100, threshold=2)
    assert "kmer_size must be larger than window_size" in str(excinfo.value)
