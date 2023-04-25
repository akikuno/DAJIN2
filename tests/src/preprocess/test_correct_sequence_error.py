# from collections import defaultdict
# from src.DAJIN2.core.preprocess.correct_sequence_error import _remove_minor_indels


# def test_remove_minor_indels():
#     count_indels = {"ins": [3, 4, 1, 2, 2], "del": [2, 2, 30, 4, 1], "sub": [1, 2, 1, 1, 1]}
#     test = _remove_minor_indels(count_indels, coverage = 1000)
#     answer = {'ins': [1, 1, 1, 1, 1], 'del': [1, 1, 30, 1, 1], 'sub': [1, 1, 1, 1, 1]}
#     assert test == answer
