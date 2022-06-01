from src.DAJIN2.preprocess import mappy_align

reffa = "tests/data/mappy/tyr_control.fa"
quefq = "tests/data/mappy/tyr_query.fq"

SAM = mappy_align.to_sam(reffa, quefq)

with open("tests/data/mappy/tyr_query.sam") as f:
    SAM_EVAL = f.read().splitlines()


def test_mappy_align():
    assert SAM == SAM_EVAL


# Create test data

# # Output as BAM
# with open("tests/data/mappy/tyr_query.sam", "w") as f:
#     f.write("\n".join(SAM))

# import pysam
# pysam.sort("-o", "tmp.bam", "tests/data/mappy/tyr_query.sam", catch_stdout=False)
# pysam.index("tmp.bam")


# END
