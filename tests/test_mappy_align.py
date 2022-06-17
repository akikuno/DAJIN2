from pathlib import Path
from src.DAJIN2.preprocess import mappy_align


def test_to_sam():
    reffa = Path("tests", "data", "mappy", "tyr_control.fa")
    quefq = Path("tests", "data", "mappy", "tyr_query.fq")
    value = mappy_align.to_sam(str(reffa), str(quefq))
    answer = Path("tests", "data", "mappy", "tyr_query.sam").read_text().strip().split("\n")
    assert value == answer


def test_remove_unmapped_reads():
    sam = Path("tests", "data", "mappy", "unmapped.sam").read_text().strip().split("\n")
    value = mappy_align.remove_unmapped_reads(sam)
    answer = Path("tests", "data", "mappy", "unmapped_removed.sam").read_text().strip().split("\n")
    assert value == answer


def test_remove_long_softclipped_reads():
    sam = Path("tests", "data", "mappy", "long_softclip.sam").read_text().strip().split("\n")
    value = mappy_align.remove_long_softclipped_reads(sam)
    answer = Path("tests", "data", "mappy", "long_softclip_removed.sam").read_text().strip().split("\n")
    assert value == answer


# Create test data

# # Output as BAM
# with open("tests/data/mappy/tyr_query.sam", "w") as f:
#     f.write("\n".join(SAM))

# import pysam
# pysam.sort("-o", "tmp.bam", "tests/data/mappy/tyr_query.sam", catch_stdout=False)
# pysam.index("tmp.bam")


# END
