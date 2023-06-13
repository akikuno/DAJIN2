# import pytest
# from src.DAJIN2.core.preprocess import correct_knockin
# from importlib import reload
# from pathlib import Path

# reload(correct_knockin)


# def test_extract_knockin_loci() -> None:
#     TEMPDIR = Path("tests/data/preprocess/correct_knockin/")
#     test = correct_knockin.extract_knockin_loci(TEMPDIR)
#     answer = {"control": {59, 60}, "random": {64, 65, 52, 53}}
#     assert test == answer


# def test_get_5mer_of_sequence() -> None:
#     sequence = "ACGTACGT"
#     test = correct_knockin.get_5mer_of_sequence(sequence)
#     answer = {2: "ACGTA", 3: "CGTAC", 4: "GTACG"}
#     assert test == answer


# def test_get_idx_of_similar_5mers() -> None:
#     kmer = {1738: "CATAA", 2451: "ATAAC"}
#     seq = {2796: "CTAAT", 2797: "TAATG"}
#     loci = {1000, 2000}
#     test = correct_knockin.get_idx_of_similar_5mers(kmer, seq, loci, n=5)
#     answer = {1738: {2796, 2797}, 2451: {2796, 2797}}
#     assert dict(test) == answer


# def test_count_indel_5mer() -> None:
#     cssplits_transposed = [
#         ["=A", "=A", "=A", "=A", "=A", "=A"],
#         ["=A", "*AC", "=A", "=A", "=A", "=A"],
#         ["=A", "-A", "=A", "=A", "=A", "=A"],
#         ["=A", "+T|+T|A", "=A", "=A", "=A", "=A"],
#         ["*AC", "*AC", "*AC", "*AC", "*AC", "*AC", "*AC"],
#         ["=A", "*AC", "-A", "+T|T|=A", "=A", "=A"],
#     ]
#     indexs = [3]
#     test = correct_knockin.count_indel_5mer(cssplits_transposed, indexs)
#     answer = {3: {"ins": [1, 1, 2, 1, 2], "del": [1, 2, 1, 1, 2], "sub": [2, 1, 1, 8, 2]}}
#     assert dict(test) == answer
