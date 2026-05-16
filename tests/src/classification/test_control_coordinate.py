from __future__ import annotations

from pathlib import Path

from DAJIN2.core.classification import control_coordinate
from DAJIN2.utils import fastx_handler


def test_rank_alleles_preserves_stx2_deletion_signature():
    signatures = {
        "control": ["=A", "=C", "=G", "=T", "=A", "=C"],
        "deletion_260": ["=A", "-C", "=G", "=T", "=A", "=C"],
        "deletion_620": ["=A", "-C", "-G", "-T", "=A", "=C"],
        "deletion_1800": ["=A", "-C", "-G", "-T", "-A", "-C"],
    }

    ranked = control_coordinate.rank_alleles("=A,-C,-G,-T,=A,=C", signatures)

    assert ranked[0][0] == "deletion_620"


def test_rank_alleles_supports_example_flox_insertion_signature():
    fasta_alleles = fastx_handler.dictionize_allele(Path("examples", "example_flox", "cables2_flox.fa"))
    assert {"control", "flox"}.issubset(fasta_alleles)

    signatures = {
        "control": ["=A", "=C", "=G", "=T"],
        "flox": ["=A", "+T|+T|+A|=C", "=G", "=T"],
    }

    ranked = control_coordinate.rank_alleles("=A,+T|+T|+A|=C,=G,=T", signatures)

    assert ranked[0][0] == "flox"


def test_rank_alleles_prefers_flox_when_insertion_payload_has_small_error():
    left_insert = "+A|+T|+A|+A|+C"
    left_insert_with_one_base_deletion = "+A|+T|+A|+C"
    right_insert = "+G|+G|+C|+T|+A"
    signatures = {
        "control": ["=A", "=C", "=G", "=T"],
        "flox": ["=A", left_insert + "|=C", "=G", right_insert + "|=T"],
        "right_insertion_only": ["=A", "=C", "=G", right_insert + "|=T"],
    }

    ranked = control_coordinate.rank_alleles(
        f"=A,{left_insert_with_one_base_deletion}|=C,=G,{right_insert}|=T",
        signatures,
    )

    assert ranked[0][0] == "flox"


def test_rank_alleles_supports_example_inversion_signature():
    fasta_alleles = fastx_handler.dictionize_allele(Path("examples", "example_inversion", "cables2_inversion.fa"))
    assert "control" in fasta_alleles
    assert Path("examples", "example_inversion", "inversion", "barcode47.fq.gz").exists()

    signatures = {
        "control": ["=A", "=C", "=G", "=T"],
        "inversion": ["=A", "=c", "=g", "=T"],
    }

    ranked = control_coordinate.rank_alleles("=A,=c,=g,=T", signatures)

    assert ranked[0][0] == "inversion"


def test_classify_records_by_control_coordinate_scores_all_alleles():
    records = [
        {"QNAME": "read1", "MIDSV": "=A,-C,-G,=T"},
        {"QNAME": "read2", "MIDSV": "=A,+G|=C,=G,=T"},
    ]
    signatures = {
        "control": ["=A", "=C", "=G", "=T"],
        "deletion": ["=A", "-C", "-G", "=T"],
        "insertion": ["=A", "+G|=C", "=G", "=T"],
    }

    classified, assignments = control_coordinate.classify_records_by_control_coordinate(
        records, signatures, no_filter=True
    )

    assert assignments == {"read1": "deletion", "read2": "insertion"}
    assert [(record["QNAME"], record["ALLELE"]) for record in classified] == [
        ("read1", "deletion"),
        ("read2", "insertion"),
    ]
