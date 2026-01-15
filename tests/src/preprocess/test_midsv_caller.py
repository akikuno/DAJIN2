from __future__ import annotations

import pytest

from DAJIN2.core.preprocess.alignment.midsv_caller import (
    convert_consecutive_indels_to_match,
    convert_flag_to_strand,
    replace_internal_n_to_d,
)

###########################################################
# replace n to d
###########################################################


def test_replace_internal_n_to_d():
    sequence = "XYZABCDEF"

    # Existing test cases
    midsv_samples = [{"MIDSV": "=N,=N,A,B,=N,=N"}]
    result = list(replace_internal_n_to_d(midsv_samples, sequence))
    assert result == [{"MIDSV": "=N,=N,A,B,=N,=N"}]

    midsv_samples = [{"MIDSV": "A,B,=N,=N,C,D"}]
    result = list(replace_internal_n_to_d(midsv_samples, sequence))
    assert result == [{"MIDSV": "A,B,-Z,-A,C,D"}]

    # New test cases
    # 1. Replace when Ns are not consecutive
    midsv_samples = [{"MIDSV": "A,=N,B,=N,C"}]
    result = list(replace_internal_n_to_d(midsv_samples, sequence))
    assert result == [{"MIDSV": "A,-Y,B,-A,C"}]

    # 2. Keep Ns when they are at the start and end
    midsv_samples = [{"MIDSV": "=N,=N,=N,A,B,C,=N,=N"}]
    result = list(replace_internal_n_to_d(midsv_samples, sequence))
    assert result == [{"MIDSV": "=N,=N,=N,A,B,C,=N,=N"}]

    # 3. Replace when Ns appear in multiple consecutive groups
    midsv_samples = [{"MIDSV": "A,B,=N,=N,C,=N,=N,D"}]
    result = list(replace_internal_n_to_d(midsv_samples, sequence))
    assert result == [{"MIDSV": "A,B,-Z,-A,C,-C,-D,D"}]

    # 4. Replace when there is a single N
    midsv_samples = [{"MIDSV": "A,B,=N,C,D"}]
    result = list(replace_internal_n_to_d(midsv_samples, sequence))
    assert result == [{"MIDSV": "A,B,-Z,C,D"}]

    # 5. Keep when there are no Ns
    midsv_samples = [{"MIDSV": "A,B,C,D,E"}]
    result = list(replace_internal_n_to_d(midsv_samples, sequence))
    assert result == [{"MIDSV": "A,B,C,D,E"}]


def test_replace_internal_n_to_d_multiple_samples():
    midsv_sample = [{"MIDSV": "=N,=N,=N,=A,=N,=C,=N,=N"}, {"MIDSV": "=N,=N,=A,=N,=N,=C,=C,=N"}]
    sequence = "GCAACCCC"
    test = replace_internal_n_to_d(midsv_sample, sequence)
    test = list(test)
    answer = [{"MIDSV": "=N,=N,=N,=A,-C,=C,=N,=N"}, {"MIDSV": "=N,=N,=A,-A,-C,=C,=C,=N"}]
    assert test == answer


def test_replace_internal_n_to_d_large_n():
    midsv_sample = [{"MIDSV": "=N,=N,=N,=N,=N,=N,=C,=N,=A"}]
    sequence = "GCAACCCCA"
    test = replace_internal_n_to_d(midsv_sample, sequence)
    test = list(test)
    answer = [{"MIDSV": "=N,=N,=N,=N,=N,=N,=C,-C,=A"}]
    assert test == answer


###########################################################
# convert_flag_to_strand
###########################################################
@pytest.mark.parametrize(
    "input_sample, expected_output",
    [
        # FLAG is 16
        ([{"FLAG": 16}], [{"STRAND": "-"}]),
        # FLAG is 2064
        ([{"FLAG": 2064}], [{"STRAND": "-"}]),
        # FLAG is neither of the above
        ([{"FLAG": 32}], [{"STRAND": "+"}]),
        # Test multiple samples
        ([{"FLAG": 16}, {"FLAG": 2064}, {"FLAG": 32}], [{"STRAND": "-"}, {"STRAND": "-"}, {"STRAND": "+"}]),
    ],
)
def test_convert_flag_to_strand(input_sample, expected_output):
    result = list(convert_flag_to_strand(iter(input_sample)))
    assert result == expected_output


###########################################################
# convert_consecutive_indels_to_match
###########################################################


@pytest.mark.parametrize(
    "cons, expected",
    [
        # simple case
        ("=A,-C,-G,-T,+T|=A", "=A,-C,-G,=T,=A"),
        ("=A,-C,-G,-T,+G|+T|=A", "=A,-C,=G,=T,=A"),
        ("=A,-C,-G,-T,+C|+G|+T|=A", "=A,=C,=G,=T,=A"),
        # no change
        ("=A,=C,=G,=T,=A", "=A,=C,=G,=T,=A"),
        # empty
        ("", ""),
    ],
)
def test_convert_consecutive_indels_to_match(cons, expected):
    assert convert_consecutive_indels_to_match(cons) == expected
