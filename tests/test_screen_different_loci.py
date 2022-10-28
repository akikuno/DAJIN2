from src.DAJIN2.core.clustering import screen_diffloci
from importlib import reload

reload(screen_diffloci)


def test_replaceNtoMatch_left():
    sample_cs = ["N", "N", "=A"]
    test = screen_diffloci.replaceNtoMatch(sample_cs)
    answer = ["=", "=", "=A"]
    assert test == answer


def test_replaceNtoMatch_right():
    sample_cs = ["=A", "N", "N"]
    test = screen_diffloci.replaceNtoMatch(sample_cs)
    answer = ["=A", "=", "="]
    assert test == answer


def test_make_table():
    cssplit_sample = [
        "*AG,=A,=A",
        "-A,=A,=A",
        "*AG,=A,=A",
    ]
    cssplit_control = [
        "=A,=A,+A|+C|=A",
        "=A,=A,=A",
        "*AG,=A,=A",
    ]
    test = screen_diffloci.make_table(cssplit_sample, cssplit_control)
    answer = [[2, 3], [1, 1], [1, 1]], [[1, 2], [1, 1], [1, 2]]
    assert test == answer


def test_screen_different_loci_repeat_insertion():
    cssplit_sample = ["+C|+C|+C|+C|+C|=A,=A,=A,=A,=C,=C,=C,=C"] * 50
    cssplit_control = ["=A,=A,=A,=A,=C,=C,=C,=C"] * 50
    sequence = "AAAACCCC"
    masks_control = [False] * len(sequence)
    alpha = 0.5
    threshold = 0.0
    test, _ = screen_diffloci.screen_different_loci(
        cssplit_sample, cssplit_control, sequence, masks_control, alpha, threshold
    )
    answer = [[0, 8.112229822253404]]
    assert test == answer


def test_screen_different_loci_no_repeat():
    cssplit_sample = ["*AG,=A,=A,=C,=C,=C"] * 100
    cssplit_control = ["=A,=A,=A,=C,=C,=C"] * 100
    sequence = "AAACCC"
    masks_control = [False] * len(sequence)
    alpha = 0.5
    threshold = 0.0
    test, _ = screen_diffloci.screen_different_loci(
        cssplit_sample, cssplit_control, sequence, masks_control, alpha, threshold
    )
    answer = [[0, 15.954589770191003]]
    assert test == answer
