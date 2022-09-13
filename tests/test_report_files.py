from src.DAJIN2.report import report_files
from importlib import reload

reload(report_files)


def test_call_fasta():
    header = "test_sequence"
    cons_seq = "A" * 81
    test = report_files.to_fasta(header, cons_seq)
    answer = ">test_sequence\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\nA\n"
    assert test == answer
