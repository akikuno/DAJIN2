import mappy as mp
import cstag
import re

reffa = "tests/data/mappy/ref.fa"
quefq = "tests/data/mappy/query.fq"

refname, refseq, _ = list(mp.fastx_read(reffa))[0]
ref = mp.Aligner(reffa)

que_map = []
for name, seq, qual in mp.fastx_read(quefq):
    for hit in ref.map(seq, cs=True):
        que_map.append([hit.r_st, hit.cigar_str, hit.cs])


def generate_query_seq(CSTAG, REFSEQ):
    cs = re.split(r"([-:+*])", CSTAG)[1:]
    cs = [i + j for i, j in zip(cs[0::2], cs[1::2])]
    query_seq = []
    start = 0
    for c in cs:
        if c[0] == ":":
            end = start + int(c[1:])
            query_seq.append(REFSEQ[start:end])
            start = end
        elif c[0] == "+":
            query_seq.append(c[1:])
        elif c[0] == "-":
            start += len(c[1:])
        elif c[0] == "*":
            query_seq.append("N")
            start += 1
    return "".join(query_seq).upper()


cslengthen = []
for r_st, cigar, cs in que_map:
    r_seq = refseq[r_st:]
    queseq = generate_query_seq(cs, r_seq)
    cslengthen.append(cstag.lengthen(cs, cigar, queseq))


###############################################################################
# Evaluate
###############################################################################

sam = "tests/data/mappy/align_cslong.sam"
with open(sam) as f:
    sam = f.read().splitlines()

cslong = [s.split("\t")[-2] for s in sam if "cs:Z" in s]


def test_cslengthen():
    for cs, cs_ans in zip(cslengthen, cslong):
        assert cs == cs_ans


# END
