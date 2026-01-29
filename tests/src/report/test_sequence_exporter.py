from __future__ import annotations

from src.DAJIN2.core import report

########################################################################
# Convert to files
########################################################################


def test_convert_to_fasta():
    header = "test_sequence"
    cons_seq = "A" * 81
    test = report.sequence_exporter.convert_to_fasta(header, cons_seq)
    answer = ">test_sequence\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\nA\n"
    assert test == answer
