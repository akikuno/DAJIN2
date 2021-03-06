from importlib import reload
from pathlib import Path
from src.DAJIN2.preprocess import mappy_align

reload(mappy_align)


def test_to_sam():
    reffa = Path("tests", "data", "mappy", "tyr_control.fa")
    quefq = Path("tests", "data", "mappy", "tyr_query.fq")
    value = mappy_align.to_sam(str(reffa), str(quefq))
    value = list(value)
    answer = Path("tests", "data", "mappy", "tyr_query.sam").read_text().strip().split("\n")
    assert value == answer


# def test_remove_unmapped_reads():
#     sam = Path("tests", "data", "mappy", "unmapped.sam").read_text().strip().split("\n")
#     value = mappy_align.remove_unmapped_reads(sam)
#     value = list(value)
#     answer = Path("tests", "data", "mappy", "unmapped_removed.sam").read_text().strip().split("\n")
#     assert value == answer


# def test_remove_overlaped_reads():
#     sam = Path("tests", "data", "mappy", "overlapped.sam").read_text().strip().split("\n")
#     value = mappy_align.remove_overlapped_reads(sam)
#     value = list(value)
#     answer = Path("tests", "data", "mappy", "overlapped_removed.sam").read_text().strip().split("\n")
#     assert value == answer


# def test_remove_overlaped_reads_2():
#     """Input is truely inversion so the function must return the same of the input.
#     """
#     sam = Path("tests", "data", "mappy", "inversion.sam").read_text().strip().split("\n")
#     value = mappy_align.remove_overlapped_reads(sam)
#     value = list(value)
#     assert value == sam


# Create test data

# # Output as BAM
# with open("tests/data/mappy/tyr_query.sam", "w") as f:
#     f.write("\n".join(SAM))

# import pysam
# pysam.sort("-o", "tmp.bam", "tests/data/mappy/tyr_query.sam", catch_stdout=False)
# pysam.index("tmp.bam")


# END
