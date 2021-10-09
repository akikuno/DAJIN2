import os
import re
import collections
from typing import NamedTuple
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


def get_reflen(sqheader: str):
    sn = [_ for _ in sqheader.split("\t") if "SN:" in _][0]
    sn = sn.replace("SN:", "")
    ln = [_ for _ in sqheader.split("\t") if "LN:" in _][0]
    ln = int(ln.replace("LN:", ""))
    return {sn: int(ln)}


def format_cstag(cstag: str):
    cstag = cstag.replace("cs:Z:", "")
    cstag = cstag.replace("-", "=D")
    cstag = cstag.replace("+", "=I")
    cstag = re.sub("[ACGT]", "M", cstag)
    cstag = re.sub("\\*[acgt][acgt]", "=S", cstag)
    cstags = cstag.split("=")
    cstags = [cs for cs in cstags if cs != '']
    return cstags


def append_next_mids_to_ins(cstags: list):
    for i, cs in enumerate(cstags):
        if "I" in cs:
            cstags[i] = cs + cstags[i+1][0]
            cstags[i+1] = cstags[i+1][1:]
    return cstags


def to_fixed_length(cs):
    if "D" in cs:
        cs = re.sub("[acgt]", "D,", cs[1:])
    elif "I" in cs:
        cs = f"{len(cs)-2}{cs[-1]},"
    elif "S" in cs:
        cs = "S,"
    elif "M" in cs:
        cs = cs.replace("M", "M,")
    return(cs)


def to_mids(cstag: str):
    cstags = format_cstag(cstag)
    cstags = append_next_mids_to_ins(cstags)
    cstags_fixlen = list(map(to_fixed_length, cstags))
    mids = ''.join(cstags_fixlen)
    return mids


def padding(mids: str, pos: int, reflen: int):
    midslen = mids.count(",")
    left_pad = "=," * (int(pos) - 1)
    right_pad = "=," * (reflen - midslen - int(pos) + 1)
    mids_padding = ''.join([left_pad, mids, right_pad])
    return mids_padding


SAMDIR = os.path.join("tests", "samTomids", "input")

SAMFILES = [os.path.join(SAMDIR, _) for _ in os.listdir(SAMDIR)]


samfile = SAMFILES[0]
# samfile = SAMFILES[6]

with open(samfile, "r") as f:
    sam = f.readlines()

sqheaders = [s for s in sam if s.startswith("@SQ")]
SQ = {}
for s in sqheaders:
    SQ.update(get_reflen(s))

contents = [s for s in sam if not s.startswith("@")]

qnames = [s.split("\t")[0] for s in contents]
qnames_counts = collections.Counter(qnames)
qnames_duplicated = [k for k, v in qnames_counts.items() if v > 1]
qnames_largedel = [k for k, v in qnames_counts.items() if v == 2]
qnames_largeinv = [k for k, v in qnames_counts.items() if v == 3]

# len(set(qnames_duplicated)) == len(set(qnames_duplicated + qnames))
for qname in qnames_largedel:
    reads = [s for s in contents if s.startswith(qname)]
    saminfo = list()
    for read in reads:
        record = read.split("\t")
        samdict = dict(
            qname=record[0],
            flag=record[1],
            reflen=int(SQ[record[2]]),
            pos=int(record[3]),
            cstag=[_ for _ in record if "cs:Z:" in _][0],
        )
        saminfo.append(samdict)
    saminfo = sorted(saminfo, key=lambda x: x['pos'])
    mids = [to_mids(s["cstag"]) for s in saminfo]
    _ = [saminfo[i].update({"mids": s}) for i, s in enumerate(mids)]
    left_len = saminfo[0]["mids"].count(",") - 1
    del_len = saminfo[1]["pos"] - saminfo[0]["pos"] - left_len
    del_seq = "D," * del_len
    mids = ''.join([saminfo[0]["mids"], del_seq, saminfo[1]["mids"]])
    mids_padding = padding(mids, saminfo[0]["pos"], saminfo[0]["reflen"])
    output = ','.join([qname, mids_padding]).rstrip(",")

with open(samfile, "r") as f:
    for line in f:
        if "@SQ" in line:
            SQ = get_reflen(line)
        if "cs:Z" in line:
            record = line.split()
            qname = record[0].replace(",", "_")
            flag = record[1]
            rname = record[2]
            pos = int(record[3])
            reflen = int(SQ[rname])
            cstag = [_ for _ in record if "cs:Z:" in _][0]
            mids = to_mids(cstag)
            mids_padding = padding(mids, pos, reflen)
            output = ','.join([qname, mids_padding]).rstrip(",")
            print(output)
            print(output.count(","))


len("MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM")
len("MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM")

"M," * 50 + "D," * 1300 + "M," * 150
