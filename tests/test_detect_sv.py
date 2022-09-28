from importlib import reload
from src.DAJIN2.core import classification

reload(classification)


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
