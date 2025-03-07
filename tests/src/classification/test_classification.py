import pytest

from DAJIN2.core.classification import classifier


@pytest.mark.parametrize(
    "input_str, expected_output",
    [
        ("=A,=C,=G,=T", 4),  # perfect match
        ("=A,=C,+T|+T|=G,=T", 2),  # insertion
        ("=A,-C,-G,=T", 0),  # deletion
        ("=A,*CT,=G,=T", 3),  # substitution
        ("=A,N,=G,=T", 3),  # unknown
    ],
)
def test_calc_match(input_str, expected_output):
    result = classifier.calc_match(input_str)
    assert result == expected_output


def test_extract_alleles_with_max_score():
    score_of_each_alleles = [
        {"QNAME": "read1", "CSSPLIT": "=A,=C,=G,=T", "SCORE": 1.0, "ALLELE": "allele1"},
        {"QNAME": "read1", "CSSPLIT": "=A,=C,=G,=T", "SCORE": 0.5, "ALLELE": "allele2"},
    ]
    result = classifier.extract_alleles_with_max_score(score_of_each_alleles)
    expected_output = [
        {"QNAME": "read1", "CSSPLIT": "=A,=C,=G,=T", "ALLELE": "allele1"},
    ]
    assert result == expected_output


def test_extract_alleles_with_max_score_multiple_reads():
    score_of_each_alleles = [
        {"QNAME": "read1", "CSSPLIT": "=A,=C,=G,=T", "SCORE": 1.0, "ALLELE": "allele1"},
        {"QNAME": "read1", "CSSPLIT": "=A,=C,=G,=T", "SCORE": 0.5, "ALLELE": "allele2"},
        {"QNAME": "read2", "CSSPLIT": "=A,=C,=G,=T", "SCORE": 0.5, "ALLELE": "allele1"},
        {"QNAME": "read2", "CSSPLIT": "=A,=C,=G,=T", "SCORE": 1.0, "ALLELE": "allele2"},
    ]
    result = classifier.extract_alleles_with_max_score(score_of_each_alleles)
    expected_output = [
        {"QNAME": "read1", "CSSPLIT": "=A,=C,=G,=T", "ALLELE": "allele1"},
        {"QNAME": "read2", "CSSPLIT": "=A,=C,=G,=T", "ALLELE": "allele2"},
    ]
    assert result == expected_output
