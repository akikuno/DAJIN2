from src.DAJIN2.clustering import module_make_scores
from importlib import reload

reload(module_make_scores)


def test_extract_cssplit_at_diffloci():
    cssplit_sample = ["-A,*AC,=A", "=A,-A,=A", "-A,N,-A"]
    diffloci = [1]
    test = module_make_scores.extract_cssplit_at_diffloci(cssplit_sample, diffloci)
    answer = [["*AC"], ["-A"], ["N"]]
    assert test == answer


def test_summation_scores():
    cssplit_sample = ["+C|+C|=A,-A,N", "=A,-A,-A", "-A,-A,*AC"]
    diffloci = [0, 1, 2]
    cssplit_diffloci = module_make_scores.extract_cssplit_at_diffloci(cssplit_sample, diffloci)
    test = module_make_scores.summation_scores(cssplit_diffloci)
    answer = [[1, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 3, 0, 0], [0, 0, 0, 0, 1, 1, 1]]
    assert test == answer


def test_repeat_del():
    cssplit_sample = ["-A,-A,=A", "=A,-A,=A", "-A,-A,-A"]
    cssplit_control = ["=A,=A,=A", "=A,=A,=A", "=A,=A,=A"]
    diffloci = [0, 1, 2]
    test = module_make_scores.make_scores(cssplit_sample, cssplit_control, diffloci)
    test = list(test)
    answer = [
        [[0, 0, 0, 0, 4, 0, 0], [0, 0, 0, 0, 6, 0, 0], [0, 0, 0, 0, 0, 0, 0]],
        [[0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 3, 0, 0], [0, 0, 0, 0, 0, 0, 0]],
        [[0, 0, 0, 0, 6, 0, 0], [0, 0, 0, 0, 9, 0, 0], [0, 0, 0, 0, 3, 0, 0]],
    ]
    assert test == answer


def test_repeat_substitution():
    cssplit_sample = ["*AC,*AC,=A", "=A,*AC,=A", "*AC,*AC,*AC"]
    cssplit_control = ["=A,=A,=A", "=A,=A,=A", "=A,=A,=A"]
    diffloci = [0, 1, 2]
    test = module_make_scores.make_scores(cssplit_sample, cssplit_control, diffloci)
    test = list(test)
    answer = [
        [[0, 0, 0, 0, 0, 4, 0], [0, 0, 0, 0, 0, 6, 0], [0, 0, 0, 0, 0, 0, 0]],
        [[0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 3, 0], [0, 0, 0, 0, 0, 0, 0]],
        [[0, 0, 0, 0, 0, 6, 0], [0, 0, 0, 0, 0, 9, 0], [0, 0, 0, 0, 0, 3, 0]],
    ]
    assert test == answer


def test_repeat_N():
    cssplit_sample = ["N,N,=A", "=A,N,=A", "N,N,N"]
    cssplit_control = ["=A,=A,=A", "=A,=A,=A", "=A,=A,=A"]
    diffloci = [0, 1, 2]
    test = module_make_scores.make_scores(cssplit_sample, cssplit_control, diffloci)
    test = list(test)
    answer = [
        [[0, 0, 0, 0, 0, 0, 4], [0, 0, 0, 0, 0, 0, 6], [0, 0, 0, 0, 0, 0, 0]],
        [[0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 3], [0, 0, 0, 0, 0, 0, 0]],
        [[0, 0, 0, 0, 0, 0, 6], [0, 0, 0, 0, 0, 0, 9], [0, 0, 0, 0, 0, 0, 3]],
    ]
    assert test == answer


def test_inversion():
    cssplit_sample = ["-a,*AC,-a", "=A,+c|c|c|=a,=A", "=A,+c|c|c|*ac,=A", "=A,*ac,=A", "=A,*ac,-a"]
    cssplit_control = ["=A,=A,=A", "=A,=A,=A", "=A,=A,=A"]
    diffloci = [0, 1, 2]
    test = module_make_scores.make_scores(cssplit_sample, cssplit_control, diffloci)
    test = list(test)
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
    cssplit_control = ["=A,=A,=A", "=A,=A,=A", "=A,=A,=A"]
    diffloci = [0, 1, 2]
    test = module_make_scores.make_scores(cssplit_sample, cssplit_control, diffloci)
    test = list(test)
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
    cssplit_control = ["=A,=A,=A", "=A,=A,=A", "=A,=A,=A"]
    diffloci = [1]
    test = module_make_scores.make_scores(cssplit_sample, cssplit_control, diffloci)
    test = list(test)
    answer = [[[0, 0, 0, 0, 3, 0, 0]], [[0, 0, 0, 0, 3, 0, 0]], [[0, 0, 0, 0, 3, 0, 0]]]
    assert test == answer


# def test_control_subtraction_equals():
#     cssplit_sample = ["-A,-A,=A", "=A,-A,=A", "-A,-A,-A"]
#     cssplit_control = ["-A,-A,=A", "=A,-A,=A", "-A,-A,-A"]
#     diffloci = [0, 1, 2]
#     test = module_make_scores.make_scores(cssplit_sample, cssplit_control, diffloci)
#     test = list(test)
#     answer = [
#         [[0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]],
#         [[0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]],
#         [[0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]],
#     ]
#     assert test == answer


# def test_control_subtraction_allzero():
#     cssplit_sample = ["-A", "-A", "-A"]
#     cssplit_control = ["-A", "-A", "-A"]
#     diffloci = [0]
#     test = module_make_scores.make_scores(cssplit_sample, cssplit_control, diffloci)
#     test = list(test)
#     answer = [[[0, 0, 0, 0, 0, 0, 0]], [[0, 0, 0, 0, 0, 0, 0]], [[0, 0, 0, 0, 0, 0, 0]]]
#     assert test == answer

