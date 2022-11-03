from src.DAJIN2.core.clustering import merge_alleles

###############################################################################
# test to replace_N_to_match
###############################################################################


def test_replace_N_to_match_left():
    cssplit = ["N", "N", "=A"]
    sequence = "GCA"
    test = merge_alleles.replace_N_to_match(cssplit, sequence)
    answer = ["=G", "=C", "=A"]
    assert test == answer


def test_replace_N_to_match_right():
    cssplit = ["=G", "N", "N"]
    sequence = "GCA"
    test = merge_alleles.replace_N_to_match(cssplit, sequence)
    answer = ["=G", "=C", "=A"]
    assert test == answer


def test_replace_N_to_match_both():
    cssplit = ["N", "N", "=A", "N"]
    sequence = "GCAT"
    test = merge_alleles.replace_N_to_match(cssplit, sequence)
    answer = ["=G", "=C", "=A", "=T"]
    assert test == answer


###############################################################################
# test to replace_maskloci_to_match
###############################################################################


def test_replace_maskloci_to_match():
    cssplit = ["=G", "D", "=A"]
    sequence = "GCA"
    maskloci = [False, True, False]
    test = merge_alleles.replace_maskloci_to_match(cssplit, sequence, maskloci)
    answer = ["=G", "=C", "=A"]
    assert test == answer
