# from __future__ import annotations
# from pathlib import Path
# from src.DAJIN2.core.clustering.past import find_knockin_loci


# def test_find_knockin_loci():
#     TEMPDIR = Path("tests/data/clustering/find_knockin_loci")
#     FASTA_ALLELES = eval(Path("tests/data/clustering/find_knockin_loci/design_cables2.jsonl").read_text())
#     CONTROL_NAME = "barcode42"
#     knockin_loci = find_knockin_loci.find_knockin_loci(TEMPDIR, FASTA_ALLELES, CONTROL_NAME)
#     test = knockin_loci["flox"]
#     answer = eval(Path("tests/data/clustering/find_knockin_loci/answer.txt").read_text())
#     assert test == answer
