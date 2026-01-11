import pytest

from DAJIN2.core.preprocess.genome_coordinate import genome_fetcher
from DAJIN2.utils.input_validator import get_available_servers


@pytest.mark.slow
def test_fetch_seq_coodinates_strand_plus():
    genome_urls = get_available_servers()

    genome = "mm39"
    gggenome_url = genome_urls["gggenome"]
    seq = "GTTAGGATTTTCAGGGTGACGACCTCCCAAGTACTCATCTGTGCAAATGT"
    test = genome_fetcher.fetch_seq_coordinates(genome, gggenome_url, seq)
    answer = {"chrom": "chr7", "strand": "+", "start": 87141775, "end": 87141825}
    assert test == answer


@pytest.mark.slow
def test_fetch_seq_coodinates_strand_minus():
    genome_urls = get_available_servers()

    genome = "mm39"
    gggenome_url = genome_urls["gggenome"]
    seq = "ACATTTGCACAGATGAGTACTTGGGAGGTCGTCACCCTGAAAATCCTAAC"
    test = genome_fetcher.fetch_seq_coordinates(genome, gggenome_url, seq)
    answer = {"chrom": "chr7", "strand": "-", "start": 87141775, "end": 87141825}
    assert test == answer


def test_fetch_seq_coordinates_multiple_regions(monkeypatch):
    genome = "hg38"
    gggenome_url = "https://gggenome.example"
    seq_subset = "GACCCTCTCTTGT"

    # fetch_bed_without_verification が複数領域を返すケース
    def fake_fetch_bed(url):
        return [
            "track header",
            ["chr1", "100", "120", "name1", "0", "+"],
            ["chr2", "200", "220", "name2", "0", "-"],
        ]

    monkeypatch.setattr(genome_fetcher, "fetch_bed_without_verification", fake_fetch_bed)

    with pytest.raises(ValueError) as ei:
        genome_fetcher.fetch_seq_coordinates(genome, gggenome_url, seq_subset)

    msg = str(ei.value)
    assert f"{gggenome_url}/{genome}/{seq_subset}" in msg
    assert "matched multiple regions" in msg
    assert "-b/--bed" in msg


def test_fetch_seq_coordinates_no_items(monkeypatch):
    genome = "hg38"
    gggenome_url = "https://gggenome.example"
    seq_subset = "GACCCTCTCTTGTGACCCTCTCTTGTGACCCTCTCTTGTGACCCTCTCTTGTGACCCTCTCTTGTGACCCTCTCTTGT"

    # fetch_bed_without_verification が "### No items found. ###" を返すケース
    def fake_fetch_bed(url):
        return [
            "track header",
            "### No items found. ###",
        ]

    monkeypatch.setattr(genome_fetcher, "fetch_bed_without_verification", fake_fetch_bed)

    with pytest.raises(ValueError) as ei:
        genome_fetcher.fetch_seq_coordinates(genome, gggenome_url, seq_subset)

    msg = str(ei.value)
    assert f"{gggenome_url}/{genome}/{seq_subset}" in msg
    assert "did not match any region" in msg
    assert "-b/--bed" in msg


@pytest.mark.slow
def test_fetch_chromosome_size():
    genome_urls = get_available_servers()
    genome = "mm39"
    chrom = "chr7"
    test = genome_fetcher.fetch_chromosome_size(genome, chrom, genome_urls["goldenpath"])
    answer = 144995196
    assert test == answer
