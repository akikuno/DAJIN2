from importlib import reload
from src.DAJIN2.clustering import find_difference
from src.DAJIN2 import clustering

reload(clustering)


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
    test = find_difference.make_table(sample_cssplit, control_cssplit)
    answer = [[1, 4], [4, 1], [4, 1]], [[3, 2], [4, 1], [3, 3]]
    assert test == answer

