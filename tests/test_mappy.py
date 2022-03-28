import re
import cstag
import mappy as mp

reffa = "tests/data/mappy/ref.fa"
quefq = "tests/data/mappy/query.fq"


refname, refseq, _ = list(mp.fastx_read(reffa))[0]


def extract_mappy_tag(REFFA, QUEFQ):
    ref = mp.Aligner(REFFA)
    qname = [qname for qname, _, _ in mp.fastx_read(QUEFQ)]
    for i, (_, seq, _) in enumerate(mp.fastx_read(QUEFQ)):
        for hit in ref.map(seq, cs=True):
            yield qname[i], hit.r_st, hit.cigar_str, hit.cs


def generate_query_seq(CSTAG, REFSEQ, R_ST):
    refseq = REFSEQ[R_ST:]
    cs = re.split(r"([-:+*])", CSTAG)[1:]
    cs = [i + j for i, j in zip(cs[0::2], cs[1::2])]
    query_seq = []
    start = 0
    for c in cs:
        if c[0] == ":":
            end = start + int(c[1:])
            query_seq.append(refseq[start:end])
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
qname_cslengthen = []
for qname, r_st, cigar, cs in extract_mappy_tag(reffa, quefq):
    queseq = generate_query_seq(cs, refseq, r_st)
    cslengthen.append(cstag.lengthen(cs, cigar, queseq))
    qname_cslengthen.append([qname, cstag.lengthen(cs, cigar, queseq)])


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


# END
