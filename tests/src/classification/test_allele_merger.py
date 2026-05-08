from DAJIN2.core.classification.allele_merger import (
    merge_minor_alleles,
    replace_negative_inf_with_most_major_allele,
)


def test_replace_negative_inf_with_most_major_allele_keeps_minor_scores():
    scores = [
        {"QNAME": "read1", "SCORE": float("-inf"), "ALLELE": "minor"},
        {"QNAME": "read1", "SCORE": -1, "ALLELE": "major"},
    ]

    result = replace_negative_inf_with_most_major_allele(scores, {"major": 10, "minor": 1})

    assert result == [
        {"QNAME": "read1", "SCORE": float("-inf"), "ALLELE": "major"},
        {"QNAME": "read1", "SCORE": -1, "ALLELE": "major"},
    ]
    assert scores[0]["ALLELE"] == "minor"


def test_merge_minor_alleles_keeps_reads_from_minor_alleles():
    scores = [
        {"QNAME": "read3", "MIDSV": "=A", "SCORE": 0, "ALLELE": "minor"},
        {"QNAME": "read1", "MIDSV": "=A", "SCORE": 0, "ALLELE": "major"},
        {"QNAME": "read2", "MIDSV": "=A", "SCORE": 0, "ALLELE": "major"},
        {"QNAME": "read3", "MIDSV": "=A", "SCORE": -1, "ALLELE": "major"},
        {"QNAME": "read1", "MIDSV": "=A", "SCORE": -1, "ALLELE": "minor"},
        {"QNAME": "read2", "MIDSV": "=A", "SCORE": -1, "ALLELE": "minor"},
    ]

    result = merge_minor_alleles(scores, threshold=2)
    read3_scores = [score for score in result if score["QNAME"] == "read3"]

    assert len(result) == len(scores)
    assert len(read3_scores) == 2
    assert {score["ALLELE"] for score in read3_scores} == {"major"}
    assert max(read3_scores, key=lambda score: score["SCORE"])["ALLELE"] == "major"
    assert scores[0]["ALLELE"] == "minor"
    assert scores[0]["SCORE"] == 0
