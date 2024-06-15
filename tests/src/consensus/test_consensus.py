from src.DAJIN2.core.consensus.consensus import (
    adjust_to_100_percent,
    call_percentage,
    replace_sequence_error,
)

###########################################################
# replace_sequence
###########################################################


def test_replace_sequence_error():
    cons_percentage = [{"A": 25, "C": 25, "SEQERROR": 50}, {"SEQERROR": 100}]
    expected_output = [{"A": 25, "C": 25}, {"N": 100}]
    assert replace_sequence_error(cons_percentage) == expected_output


###########################################################
# adjust_to_100_percent
###########################################################


def test_adjust_to_100_percent():
    test = [{"A": 25, "C": 25}, {"N": 100}]
    expected = [{"A": 50, "C": 50}, {"N": 100}]
    assert adjust_to_100_percent(test) == expected


def test_adjust_to_100_percent_float():
    test = [{"A": 20.1, "C": 19.9}]
    expected = [{"A": 50.25, "C": 49.75}]
    assert adjust_to_100_percent(test) == expected


###########################################################
# call_percentage
###########################################################


def test_call_percentage():
    cssplits = [["+A|=C", "-T", "=C", "=A", "=T"], ["-C", "=T", "=C", "*AT", "*AT"]]
    mutation_loci = [{"+", "-"}, {"-"}, {}, {}, {"*"}]
    expected_output = [
        {"+A|=C": 50.0, "-C": 50.0},
        {"-T": 50.0, "=T": 50.0},
        {"=C": 100.0},
        {"=A": 100.0},
        {"=T": 50.0, "*AT": 50.0},
    ]
    assert call_percentage(cssplits, mutation_loci) == expected_output
