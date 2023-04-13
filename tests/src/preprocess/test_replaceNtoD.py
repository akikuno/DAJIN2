from DAJIN2.core.preprocess.replace_NtoD import _replaceNtoD


def test_replaceNtoD():
    cssplits_sample = [
        ["N", "N", "N", "=A", "N", "=C", "N", "N"],
        ["N", "N", "=A", "N", "N", "=C", "=C", "N"],
    ]
    sequence = "GCAACCCC"
    test = _replaceNtoD(cssplits_sample, sequence)
    answer = [["N", "N", "N", "=A", "-C", "=C", "N", "N"], ["N", "N", "=A", "-A", "-C", "=C", "=C", "N"]]
    assert test == answer


def test_replaceNtoD_large_N():
    cssplits_sample = [
        ["N", "N", "N", "N", "N", "N", "=C", "N", "=A"],
    ]
    sequence = "GCAACCCCA"
    test = _replaceNtoD(cssplits_sample, sequence)
    answer = [["N", "N", "N", "N", "N", "N", "=C", "-C", "=A"]]
    assert test == answer
