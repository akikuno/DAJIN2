from src.DAJIN2.consensus import module_consensus
from importlib import reload

reload(module_consensus)


def test_call_fasta():
    header = "test_sequence"
    cons_percentage_by_key = [{"=A": 1.0} for _ in range(81)]
    test = module_consensus.call_fasta(header, cons_percentage_by_key)
    answer = ">test_sequence\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\nA\n"
    assert test == answer


def test_call_fasta_seq_match():
    cons_percentage_by_key = [{"=A": 1.0}, {"=T": 0.9, "-T": 0.1}]
    test = module_consensus.call_fasta_seq(cons_percentage_by_key)
    answer = "AT"
    assert test == answer


def test_call_fasta_seq_deletion():
    cons_percentage_by_key = [{"=A": 1.0}, {"-A": 0.9, "=A": 0.1}, {"=T": 1.0}]
    test = module_consensus.call_fasta_seq(cons_percentage_by_key)
    answer = "AT"
    assert test == answer


def test_call_fasta_seq_substitution():
    cons_percentage_by_key = [{"=A": 1.0}, {"*AC": 0.9, "-A": 0.1}, {"=T": 1.0}]
    test = module_consensus.call_fasta_seq(cons_percentage_by_key)
    answer = "ACT"
    assert test == answer


def test_call_fasta_seq_inversion():
    cons_percentage_by_key = [{"=A": 1.0}, {"=a": 0.9, "-a": 0.1}, {"=T": 1.0}]
    test = module_consensus.call_fasta_seq(cons_percentage_by_key)
    answer = "AaT"
    assert test == answer


def test_call_fasta_seq_insertion_match():
    cons_percentage_by_key = [{"=A": 1.0}, {"+G|+G|+G|=A": 0.9, "-A": 0.1}, {"=T": 1.0}]
    test = module_consensus.call_fasta_seq(cons_percentage_by_key)
    answer = "AGGGAT"
    assert test == answer


def test_call_fasta_seq_insertion_deletion():
    cons_percentage_by_key = [{"=A": 1.0}, {"+G|+G|+G|-A": 0.9, "-A": 0.1}, {"=T": 1.0}]
    test = module_consensus.call_fasta_seq(cons_percentage_by_key)
    answer = "AGGGT"
    assert test == answer


def test_call_fasta_seq_insertion_substitution():
    cons_percentage_by_key = [{"=A": 1.0}, {"+G|+G|+G|*AT": 0.9, "-A": 0.1}, {"=T": 1.0}]
    test = module_consensus.call_fasta_seq(cons_percentage_by_key)
    answer = "AGGGTT"
    assert test == answer


def test_call_fasta_seq_insertion_N():
    cons_percentage_by_key = [{"=A": 1.0}, {"+G|+G|+G|N": 0.9, "-A": 0.1}, {"=T": 1.0}]
    test = module_consensus.call_fasta_seq(cons_percentage_by_key)
    answer = "AGGGNT"
    assert test == answer


# def test_module_consensus():
#     cssplits = ["=A,*GA,=C", "=A,*GA,-C", "=A,*GA,-C"]
#     test = module_consensus.call(cssplits)
#     answer = "=A,*GA,-C"
#     assert test == answer

