from src.DAJIN2.core.clustering.make_scores import extract_cssplit_at_diffloci
from src.DAJIN2.core.clustering.make_scores import sum_scores
from src.DAJIN2.core.clustering.make_scores import make_scores


def test_extract_cssplit_at_diffloci():
    cssplit_sample = ["=A,*AC,=A", "=A,-A,=A", "=A,N,-A"]
    diffloci = [1]
    test = extract_cssplit_at_diffloci(cssplit_sample, diffloci)
    answer = [["*AC"], ["-A"], ["N"]]
    assert test == answer


def test_extract_cssplit_at_diffloci_2():
    cssplit_sample = ["=A,*AC,*AG", "=A,-A,*AT", "=A,N,-A"]
    diffloci = [1, 2]
    test = extract_cssplit_at_diffloci(cssplit_sample, diffloci)
    answer = [["*AC", "*AG"], ["-A", "*AT"], ["N", "-A"]]
    assert test == answer


def test_sum_scores():
    cssplit_sample = ["+C|+C|=A,-A,N", "=A,-A,-A", "-A,-A,*AC"]
    diffloci = [0, 1, 2]
    cssplit_diffloci = extract_cssplit_at_diffloci(cssplit_sample, diffloci)
    test = sum_scores(cssplit_diffloci)
    answer = [[1, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 3, 0, 0], [0, 0, 0, 0, 1, 1, 1]]
    assert test == answer


def test_repeat_del():
    cssplit_sample = ["-A,-A,=A", "=A,-A,=A", "-A,-A,-A"]
    diffloci = [0, 1, 2]
    test = make_scores(cssplit_sample, diffloci)
    answer = [
        [[0, 0, 0, 0, 4, 0, 0], [0, 0, 0, 0, 6, 0, 0], [0, 0, 0, 0, 0, 0, 0]],
        [[0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 3, 0, 0], [0, 0, 0, 0, 0, 0, 0]],
        [[0, 0, 0, 0, 6, 0, 0], [0, 0, 0, 0, 9, 0, 0], [0, 0, 0, 0, 3, 0, 0]],
    ]
    assert test == answer


def test_repeat_substitution():
    cssplit_sample = ["*AC,*AC,=A", "=A,*AC,=A", "*AC,*AC,*AC"]
    diffloci = [0, 1, 2]
    test = make_scores(cssplit_sample, diffloci)
    answer = [
        [[0, 0, 0, 0, 0, 4, 0], [0, 0, 0, 0, 0, 6, 0], [0, 0, 0, 0, 0, 0, 0]],
        [[0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 3, 0], [0, 0, 0, 0, 0, 0, 0]],
        [[0, 0, 0, 0, 0, 6, 0], [0, 0, 0, 0, 0, 9, 0], [0, 0, 0, 0, 0, 3, 0]],
    ]
    assert test == answer


def test_repeat_N():
    cssplit_sample = ["N,N,=A", "=A,N,=A", "N,N,N"]
    diffloci = [0, 1, 2]
    test = make_scores(cssplit_sample, diffloci)
    answer = [
        [[0, 0, 0, 0, 0, 0, 4], [0, 0, 0, 0, 0, 0, 6], [0, 0, 0, 0, 0, 0, 0]],
        [[0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 3], [0, 0, 0, 0, 0, 0, 0]],
        [[0, 0, 0, 0, 0, 0, 6], [0, 0, 0, 0, 0, 0, 9], [0, 0, 0, 0, 0, 0, 3]],
    ]
    assert test == answer


def test_inversion():
    cssplit_sample = ["-a,*AC,-a", "=A,+c|c|c|=a,=A", "=A,+c|c|c|*ac,=A", "=A,*ac,=A", "=A,*ac,-a"]
    diffloci = [0, 1, 2]
    test = make_scores(cssplit_sample, diffloci)
    answer = [
        [[0, 0, 0, 0, -1, 0, 0], [0, 0, 0, 0, 0, 3, 0], [0, 0, 0, 0, -2, 0, 0]],
        [[0, 0, 0, 0, 0, 0, 0], [-4, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]],
        [[0, 0, 0, 0, 0, 0, 0], [0, 0, -4, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]],
        [[0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, -3, 0], [0, 0, 0, 0, 0, 0, 0]],
        [[0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, -3, 0], [0, 0, 0, 0, -2, 0, 0]],
    ]
    assert test == answer


def test_insertion():
    cssplit_sample = ["=A,+C|C|C|=A,=A", "=A,+C|C|C|-A,=A", "=A,+C|C|C|*CG,=A", "=A,+C|C|C|N,=A", "=A,+c|c|c|n,=A"]
    diffloci = [0, 1, 2]
    test = make_scores(cssplit_sample, diffloci)
    answer = [
        [[0, 0, 0, 0, 0, 0, 0], [4, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]],
        [[0, 0, 0, 0, 0, 0, 0], [0, 4, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]],
        [[0, 0, 0, 0, 0, 0, 0], [0, 0, 4, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]],
        [[0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 4, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]],
        [[0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -4, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]],
    ]
    assert test == answer


def test_small_diffloci():
    cssplit_sample = ["-A,-A,=A", "=A,-A,=A", "-A,-A,-A"]
    diffloci = [1]
    test = make_scores(cssplit_sample, diffloci)
    answer = [[[0, 0, 0, 0, 3, 0, 0]], [[0, 0, 0, 0, 3, 0, 0]], [[0, 0, 0, 0, 3, 0, 0]]]
    assert test == answer
