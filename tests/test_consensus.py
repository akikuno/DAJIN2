from src.DAJIN2.core.consensus import consensus_main as consensus
from importlib import reload

reload(consensus)


def test_call_sequence_match():
    cons_percentage_by_key = [{"=A": 1.0}, {"=T": 0.9, "-T": 0.1}]
    test = consensus.call_sequence(cons_percentage_by_key)
    answer = "AT"
    assert test == answer


def test_call_sequence_deletion():
    cons_percentage_by_key = [{"=A": 1.0}, {"-A": 0.9, "=A": 0.1}, {"=T": 1.0}]
    test = consensus.call_sequence(cons_percentage_by_key)
    answer = "AT"
    assert test == answer


def test_call_sequence_substitution():
    cons_percentage_by_key = [{"=A": 1.0}, {"*AC": 0.9, "-A": 0.1}, {"=T": 1.0}]
    test = consensus.call_sequence(cons_percentage_by_key)
    answer = "ACT"
    assert test == answer


def test_call_sequence_inversion():
    cons_percentage_by_key = [{"=A": 1.0}, {"=a": 0.9, "-a": 0.1}, {"=T": 1.0}]
    test = consensus.call_sequence(cons_percentage_by_key)
    answer = "AaT"
    assert test == answer


def test_call_sequence_insertion_match():
    cons_percentage_by_key = [{"=A": 1.0}, {"+G|+G|+G|=A": 0.9, "-A": 0.1}, {"=T": 1.0}]
    test = consensus.call_sequence(cons_percentage_by_key)
    answer = "AGGGAT"
    assert test == answer


def test_call_sequence_insertion_deletion():
    cons_percentage_by_key = [{"=A": 1.0}, {"+G|+G|+G|-A": 0.9, "-A": 0.1}, {"=T": 1.0}]
    test = consensus.call_sequence(cons_percentage_by_key)
    answer = "AGGGT"
    assert test == answer


def test_call_sequence_insertion_substitution():
    cons_percentage_by_key = [{"=A": 1.0}, {"+G|+G|+G|*AT": 0.9, "-A": 0.1}, {"=T": 1.0}]
    test = consensus.call_sequence(cons_percentage_by_key)
    answer = "AGGGTT"
    assert test == answer


def test_call_sequence_insertion_N():
    cons_percentage_by_key = [{"=A": 1.0}, {"+G|+G|+G|N": 0.9, "-A": 0.1}, {"=T": 1.0}]
    test = consensus.call_sequence(cons_percentage_by_key)
    answer = "AGGGNT"
    assert test == answer

