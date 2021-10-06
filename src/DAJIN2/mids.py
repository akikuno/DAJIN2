import os
import re
# import dataclasses


# @dataclasses.dataclass
# class Args:
#     sample = str,
#     allele = str,
#     control = str,
#     threads = int = 1,


# args = Args()

# args.sample = "tests/data/query.fq"
# args.control = "tests/data/query.fq"
SAMDIR = os.path.join(".tmpDAJIN", "sam")

SAMFILES = [os.path.join(SAMDIR, _) for _ in os.listdir(SAMDIR)]


sam = SAMFILES[0]


with open(sam, "r") as f:
    SQ = {}
    for line in f:
        if "@SQ" in line:
            sn = [_ for _ in line.split() if "SN:" in _][0]
            sn = sn.replace("SN:", "")
            ln = [_ for _ in line.split() if "LN:" in _][0]
            ln = int(ln.replace("LN:", ""))
            SQ.update({sn: ln})
            print(SQ)
        if "cs:Z" in line:
            record = line.split()
            qname = record[0]
            flag = record[1]
            rname = record[2]
            pos = record[3]
            reflen = SQ[rname]
            cs = [_ for _ in record if "cs:Z:" in _][0]
            print(cs, reflen)
