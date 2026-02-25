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


def test_export_reference_to_fasta_keeps_control_header(tmp_path):
    tempdir = tmp_path
    sample_name = "sample"
    path_input_dir = tempdir / sample_name / "fasta"
    path_output_dir = tempdir / "report" / "FASTA" / sample_name
    path_input_dir.mkdir(parents=True)
    path_output_dir.mkdir(parents=True)

    (path_input_dir / "control.fasta").write_text(">control\nACGT\n", encoding="utf-8")
    (path_input_dir / "deletion.fasta").write_text(">deletion\nTGCA\n", encoding="utf-8")

    report.sequence_exporter.export_reference_to_fasta(tempdir, sample_name)

    assert (path_output_dir / "control.fasta").read_text(encoding="utf-8").startswith(">control\n")
    assert (path_output_dir / "deletion.fasta").read_text(encoding="utf-8").startswith(">sample_deletion\n")
