from collections import defaultdict
from src.DAJIN2.core.preprocess.replace_N_to_D import replaceNtoD


def test_replaceNtoD():
    cssplits_sample = [
        ["N", "N", "N", "=A", "N", "=C", "N", "N"],
        ["N", "N", "=A", "N", "N", "=C", "=C", "N"],
    ]
    sequence = "GCAACCCC"
    test = replaceNtoD(cssplits_sample, sequence)
    answer = [["N", "N", "N", "=A", "-C", "=C", "N", "N"], ["N", "N", "=A", "-A", "-C", "=C", "=C", "N"]]
    assert test == answer
