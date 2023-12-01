import pytest

from DAJIN2.core.preprocess.mutation_extractor import merge_surrounding_index


@pytest.mark.parametrize(
    "input_list, expected_output",
    [
        ([5, 10], {5, 6, 7, 8, 9, 10}),
        ([5, 11], {5, 11}),
        ([1, 2, 3, 10, 11, 15], {1, 2, 3, 10, 11, 12, 13, 14, 15}),
        ([0, 5, 6, 7], {0, 1, 2, 3, 4, 5, 6, 7}),
        ([2, 8, 9, 14, 19, 20, 25], {2, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25}),
        ([1, 7, 12, 18], {1, 7, 8, 9, 10, 11, 12, 18}),
        ([10, 15, 20, 26], {10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 26}),
    ],
)
def test_merge_surrounding_index(input_list, expected_output):
    assert merge_surrounding_index(input_list) == expected_output
