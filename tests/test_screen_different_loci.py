from src.DAJIN2.clustering import module_screen_diffloci
from importlib import reload

reload(module_screen_diffloci)


def test_replaceN_left():
    sample_cs = ["N", "N", "=A"]
    test = module_screen_diffloci.replaceN(sample_cs)
    answer = ["@", "@", "=A"]
    assert test == answer


def test_replaceN_right():
    sample_cs = ["=A", "N", "N"]
    test = module_screen_diffloci.replaceN(sample_cs)
    answer = ["=A", "@", "@"]
    assert test == answer


def test_make_table():
    sample_cssplit = [
        "*AG,=A,=A",
        "-A,=A,=A",
        "*AG,=A,=A",
    ]
    control_cssplit = [
        "=A,=A,+A|+C|=A",
        "=A,=A,=A",
        "*AG,=A,=A",
    ]
    test = module_screen_diffloci.make_table(sample_cssplit, control_cssplit)
    answer = [[1, 4], [4, 1], [4, 1]], [[3, 2], [4, 1], [3, 3]]
    assert test == answer


def test_screen_different_loci_repeat():
    sample_cssplit = [
        "+C|+C|+C|+C|+C|=A,=A,=A,=A,=C,=C,=C,=C",
        "+C|+C|+C|+C|+C|=A,=A,=A,=A,=C,=C,=C,=C",
        "+C|+C|+C|+C|+C|=A,=A,=A,=A,=C,=C,=C,=C",
    ]
    control_cssplit = [
        "=A,=A,=A,=A,=C,=C,=C,=C",
        "=A,=A,=A,=A,=C,=C,=C,=C",
        "=A,=A,=A,=A,=C,=C,=C,=C",
    ]
    sequence = "AAAACCCC"
    alpha = 0.1
    threshold = 0.0
    test = module_screen_diffloci.screen_different_loci(sample_cssplit, control_cssplit, sequence, alpha, threshold)
    answer = [0]
    assert test == answer


def test_screen_different_loci_no_repeat():
    sample_cssplit = [
        "*AG,=A,=A,=C,=C,=C",
        "-A,=A,=A,=C,=C,=C",
        "*AG,=A,=A,=C,=C,=C",
    ]
    control_cssplit = [
        "=A,=A,=A,=C,=C,=C",
        "=A,=A,=A,=C,=C,=C",
        "=A,=A,=A,=C,=C,=C",
    ]
    sequence = "AAACCC"
    alpha = 0.25
    threshold = 0.0
    test = module_screen_diffloci.screen_different_loci(sample_cssplit, control_cssplit, sequence, alpha, threshold)
    answer = [0]
    assert test == answer
