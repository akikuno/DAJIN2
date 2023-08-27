import pytest

from pathlib import Path
from src.DAJIN2.core.preprocess import genome_fetcher


@pytest.mark.slow
def test_fetch_seq_coodinates_strand_plus():
    genome = "mm39"
    blat_url = "https://genome.ucsc.edu/cgi-bin/hgBlat"
    seq = "GTTAGGATTTTCAGGGTGACGACCTCCCAAGTACTCATCTGTGCAAATGT"
    test = genome_fetcher._fetch_seq_coordinates(genome, blat_url, seq)
    answer = {"chrom": "chr7", "strand": "+", "start": 87141776, "end": 87141825}
    assert test == answer


@pytest.mark.slow
def test_fetch_seq_coodinates_strand_minus():
    genome = "mm39"
    blat_url = "https://genome.ucsc.edu/cgi-bin/hgBlat"
    seq = "ACATTTGCACAGATGAGTACTTGGGAGGTCGTCACCCTGAAAATCCTAAC"
    test = genome_fetcher._fetch_seq_coordinates(genome, blat_url, seq)
    answer = {"chrom": "chr7", "strand": "-", "start": 87141776, "end": 87141825}
    assert test == answer


@pytest.mark.slow
def test_fetch_seq_coodinates_error():
    genome = "mm39"
    blat_url = "https://genome.ucsc.edu/cgi-bin/hgBlat"
    seq = "XXXXXXXXXXXXXXXXX"
    with pytest.raises(ValueError) as e:
        genome_fetcher._fetch_seq_coordinates(genome, blat_url, seq)
    assert str(e.value) == f"{seq} is not found in {genome}"


@pytest.mark.slow
def test_fetch_coordinates_strand_plus():
    genome_coordinates = {"genome": "mm39"}
    genome_urls = {
        "blat": "https://genome.ucsc.edu/cgi-bin/hgBlat",
        "goldenpath": "http://hgdownload.soe.ucsc.edu/goldenPath",
    }
    seq = Path("tests", "data", "preprocess", "genome_fetcher", "sequence_plus.txt").read_text().strip()
    test = genome_fetcher.fetch_coordinates(genome_coordinates, genome_urls, seq)
    answer = {"genome": "mm39", "chrom": "chr7", "start": 87140387, "end": 87143847, "strand": "+"}
    assert test == answer


@pytest.mark.slow
def test_fetch_coordinates_strand_minus():
    genome_coordinates = {"genome": "mm39"}
    genome_urls = {
        "blat": "https://genome.ucsc.edu/cgi-bin/hgBlat",
        "goldenpath": "http://hgdownload.soe.ucsc.edu/goldenPath",
    }
    seq = Path("tests", "data", "preprocess", "genome_fetcher", "sequence_minus.txt").read_text().strip()
    test = genome_fetcher.fetch_coordinates(genome_coordinates, genome_urls, seq)
    answer = {"genome": "mm39", "chrom": "chr7", "start": 87140387, "end": 87143847, "strand": "-"}
    assert test == answer


@pytest.mark.slow
def test_fetch_chromosome_size():
    genome_coordinates = {"genome": "mm39", "chrom": "chr7"}
    genome_urls = {
        "blat": "https://genome.ucsc.edu/cgi-bin/hgBlat",
        "goldenpath": "http://hgdownload.soe.ucsc.edu/goldenPath",
    }
    test = genome_fetcher.fetch_chromosome_size(genome_coordinates, genome_urls)
    answer = 144995196
    assert test == answer
