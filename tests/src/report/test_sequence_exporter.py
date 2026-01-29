from __future__ import annotations

from pathlib import Path

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


def test_export_reference_to_fasta_uses_display_names(tmp_path):
    sample_name = "sample1"
    fasta_dir = Path(tmp_path, sample_name, "fasta")
    fasta_dir.mkdir(parents=True, exist_ok=True)

    Path(fasta_dir, "deletion001__uuid.fasta").write_text(">deletion001__uuid\nAC\n")
    Path(fasta_dir, "control.fasta").write_text(">control\nTT\n")

    sv_name_map = {"deletion001__uuid": "unintended_deletion_1"}

    report.sequence_exporter.export_reference_to_fasta(tmp_path, sample_name, sv_name_map)

    output_dir = Path(tmp_path, "report", "FASTA", sample_name)
    assert (
        Path(output_dir, "unintended_deletion_1.fasta").read_text().startswith(f">{sample_name}_unintended_deletion_1")
    )
    assert Path(output_dir, "control.fasta").read_text().startswith(f">{sample_name}_control")


def test_convert_to_html_resolves_internal_sv_name(tmp_path):
    sample_name = "sample1"
    internal_name = "deletion001__uuid"
    display_name = "unintended_deletion_1"

    midsv_dir = Path(tmp_path, sample_name, "midsv")
    midsv_dir.mkdir(parents=True, exist_ok=True)
    Path(midsv_dir, f"consensus_{internal_name}.jsonl").write_text('"=A"\n')

    sv_name_map = {internal_name: display_name}
    fasta_alleles = {internal_name: "A"}
    header = f"allele01_{display_name}_50%"
    cons_midsv_tag = ["=A"]

    html = report.sequence_exporter.convert_to_html(
        tmp_path, sample_name, fasta_alleles, header, cons_midsv_tag, sv_name_map
    )

    assert "unintended deletion 1" in html
