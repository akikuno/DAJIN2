from src.DAJIN2 import clustering
from importlib import reload

reload(clustering)


def test_repeat_del():
    sample_cssplit = ["-A,-A,=A", "=A,-A,=A", "-A,-A,-A"]
    diff_loci = [0, 1, 2]
    test = clustering.make_scores(sample_cssplit, diff_loci)
    test = list(test)
    answer = [
        [[0, 0, 0, 0, 2, 0, 0], [0, 0, 0, 0, 2, 0, 0], [0, 0, 0, 0, 0, 0, 0]],
        [[0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0]],
        [[0, 0, 0, 0, 3, 0, 0], [0, 0, 0, 0, 3, 0, 0], [0, 0, 0, 0, 3, 0, 0]],
    ]
    assert test == answer


def test_repeat_substitution():
    sample_cssplit = ["*AC,*AC,=A", "=A,*AC,=A", "*AC,*AC,*AC"]
    diff_loci = [0, 1, 2]
    test = clustering.make_scores(sample_cssplit, diff_loci)
    test = list(test)
    answer = [
        [[0, 0, 0, 0, 0, 2, 0], [0, 0, 0, 0, 0, 2, 0], [0, 0, 0, 0, 0, 0, 0]],
        [[0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0]],
        [[0, 0, 0, 0, 0, 3, 0], [0, 0, 0, 0, 0, 3, 0], [0, 0, 0, 0, 0, 3, 0]],
    ]
    assert test == answer


def test_repeat_N():
    sample_cssplit = ["N,N,=A", "=A,N,=A", "N,N,N"]
    diff_loci = [0, 1, 2]
    test = clustering.make_scores(sample_cssplit, diff_loci)
    test = list(test)
    answer = [
        [[0, 0, 0, 0, 0, 0, 2], [0, 0, 0, 0, 0, 0, 2], [0, 0, 0, 0, 0, 0, 0]],
        [[0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0]],
        [[0, 0, 0, 0, 0, 0, 3], [0, 0, 0, 0, 0, 0, 3], [0, 0, 0, 0, 0, 0, 3]],
    ]
    assert test == answer


def test_inversion():
    sample_cssplit = ["-a,*AC,-a", "=A,+c|c|c|=a,=A", "=A,+c|c|c|*ac,=A", "=A,*ac,=A", "=A,*ac,-a"]
    diff_loci = [0, 1, 2]
    test = clustering.make_scores(sample_cssplit, diff_loci)
    test = list(test)
    answer = [
        [[0, 0, 0, 0, -1, 0, 0], [0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, -1, 0, 0]],
        [[0, 0, 0, 0, 0, 0, 0], [-4, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]],
        [[0, 0, 0, 0, 0, 0, 0], [0, 0, -4, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]],
        [[0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, 0, 0]],
        [[0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, -1, 0], [0, 0, 0, 0, -1, 0, 0]],
    ]
    assert test == answer


def test_insertion():
    sample_cssplit = ["=A,+C|C|C|=A,=A", "=A,+C|C|C|-A,=A", "=A,+C|C|C|*CG,=A", "=A,+C|C|C|N,=A", "=A,+c|c|c|n,=A"]
    diff_loci = [0, 1, 2]
    test = clustering.make_scores(sample_cssplit, diff_loci)
    test = list(test)
    answer = [
        [[0, 0, 0, 0, 0, 0, 0], [4, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]],
        [[0, 0, 0, 0, 0, 0, 0], [0, 4, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]],
        [[0, 0, 0, 0, 0, 0, 0], [0, 0, 4, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]],
        [[0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 4, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]],
        [[0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -4, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]],
    ]
    assert test == answer


def test_small_diff_loci():
    sample_cssplit = ["-A,-A,=A", "=A,-A,=A", "-A,-A,-A"]
    diff_loci = [1]
    test = clustering.make_scores(sample_cssplit, diff_loci)
    test = list(test)
    answer = [[[0, 0, 0, 0, 1, 0, 0]], [[0, 0, 0, 0, 1, 0, 0]], [[0, 0, 0, 0, 1, 0, 0]]]
    assert test == answer

