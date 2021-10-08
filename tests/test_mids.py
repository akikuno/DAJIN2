import os
from typing import List, Dict
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


def get_reflen(line: str) -> Dict[str, int]:
    try:
        if SQ:
            pass
    except NameError:
        SQ = dict()
    sn = [_ for _ in line.split() if "SN:" in _][0]
    sn = sn.replace("SN:", "")
    ln = [_ for _ in line.split() if "LN:" in _][0]
    ln = int(ln.replace("LN:", ""))
    SQ.update({sn: ln})
    return SQ


def format_cstag(cstag: str) -> List[str]:
    cstag = cstag.replace("cs:Z:", "")
    cstag = cstag.replace("-", "=D")
    cstag = cstag.replace("+", "=I")
    cstag = re.sub("[ACGT]", "M", cstag)
    cstag = re.sub("\\*[acgt][acgt]", "=S", cstag)
    cstags = cstag.split("=")
    cstags = [cs for cs in cstags if cs != '']
    return cstags


def to_fixed_length(cs):
    if "D" in cs:
        cs = re.sub("[acgt]", "D,", cs[1:])
    elif "I" in cs:
        cs = f"{len(cs)-1}I,"
    elif "S" in cs:
        cs = "S,"
    elif "M" in cs:
        cs = cs.replace("M", "M,")
    return(cs)


SAMDIR = os.path.join("tests", "samTomids", "input")

SAMFILES = [os.path.join(SAMDIR, _) for _ in os.listdir(SAMDIR)]


sam = SAMFILES[3]


with open(sam, "r") as f:
    for line in f:
        if "@SQ" in line:
            SQ = get_reflen(line)
        if "cs:Z" in line:
            record = line.split()
            qname = record[0]
            flag = record[1]
            rname = record[2]
            pos = record[3]
            reflen = SQ[rname]
            cstag = [_ for _ in record if "cs:Z:" in _][0]
            cstags = format_cstag(cstag)
            cstags_fixlen = list(map(to_fixed_length, cstags))
            mids = ''.join(cstags_fixlen)
            #! PADDING =========================
            midslen = mids.count(",")
            left_pad = "=," * (int(pos) - 1)
            right_pad = "=," * (reflen - midslen - int(pos) + 1)
            mids = ''.join([left_pad, mids, right_pad])
            output = ','.join([qname, mids]).rstrip(",")
            print(output)
