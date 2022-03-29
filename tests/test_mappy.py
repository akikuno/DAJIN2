# https://github.com/lh3/minimap2/blob/master/python/README.rst

import cstag
import pysam
from src.DAJIN2 import mappy2sam

reffa = "tests/data/mappy/ref.fa"
quefq = "tests/data/mappy/query.fq"


SAM = mappy2sam(reffa, quefq)

# ------------------------------------------------------------------------------
# cslengthen
# ------------------------------------------------------------------------------

cslengthen = []
qname_cslengthen = []
for sam in SAM:
    if sam[0] == "@":
        continue
    sam = sam.split("\t")
    qname, cigar, qseq, cs = sam[0], sam[5], sam[9], sam[-1]
    cslong = cstag.lengthen(cs, cigar, qseq)
    cslengthen.append(cslong)
    qname_cslengthen.append([qname, cslong])

###############################################################################
# Evaluate
###############################################################################

sam = "tests/data/mappy/align_cslong.sam"
with open(sam) as f:
    sam = f.read().splitlines()

cslong = [s.split("\t")[-2] for s in sam if "cs:Z" in s]
qname = [q.split("\t")[0] for q in sam if "cs:Z" in q]


def test_cslengthen():
    for cs, cs_ans in zip(cslengthen, cslong):
        assert cs == cs_ans


def test_qname_cslengthen():
    for cs, cs_ans in zip(qname_cslengthen, zip(qname, cslong)):
        assert tuple(cs) == cs_ans


sam = "tests/data/mappy/mappy2sam.sam"
with open(sam) as f:
    sam = f.read().splitlines()

# pysam.sort("-o", "tests/data/mappy2sam.bam", "tmp.sam", catch_stdout=False)
# pysam.index("tests/data/mappy2sam.bam")


def test_mappy2sam():
    assert SAM == sam


# END
