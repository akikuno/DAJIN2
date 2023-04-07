from collections import defaultdict
from src.DAJIN2.core.preprocess.correct_sequence_error import remove_minor_indels


def test_remove_minor_indels():
    cssplits = [["A"]] * 1000
    count_5mer = [
        {"ins": [3, 4, 1, 2, 2], "del": [2, 2, 30, 4, 1], "sub": [1, 2, 1, 1, 1]},
        {"ins": [3, 4, 1, 2, 2], "del": [2, 2, 3, 4, 1], "sub": [1, 2, 1, 1, 1]},
    ]
    expected_result = [
        defaultdict(list, {"ins": [1, 1, 1, 1, 1], "del": [2, 2, 30, 4, 1], "sub": [1, 1, 1, 1, 1]}),
        defaultdict(list, {"ins": [1, 1, 1, 1, 1], "del": [1, 1, 1, 1, 1], "sub": [1, 1, 1, 1, 1]}),
    ]
    assert remove_minor_indels(cssplits, count_5mer) == expected_result
