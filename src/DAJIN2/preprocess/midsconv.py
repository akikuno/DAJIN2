from itertools import groupby
import re
from concurrent.futures import ProcessPoolExecutor


def get_sqheaders(sqheaders: list) -> dict:
    sq = {}
    for sqheader in sqheaders:
        sn = [_ for _ in sqheader.split("\t") if "SN:" in _][0]
        sn = sn.replace("SN:", "")
        ln = [_ for _ in sqheader.split("\t") if "LN:" in _][0]
        ln = int(ln.replace("LN:", ""))
        sq.update({sn: ln})
    return sq


def format_cstag(cstag: str) -> list:
    cstag = cstag.replace("cs:Z:", "")
    cstag = cstag.replace("-", "=D")
    cstag = cstag.replace("+", "=I")
    cstag = re.sub("[ACGT]", "M", cstag)
    cstag = re.sub("\\*[acgt][acgt]", "=S", cstag)
    cstags = cstag.split("=")
    cstags = [cs for cs in cstags if cs != ""]
    return cstags


def append_next_mids_to_ins(cstags: list) -> list:
    for i, cs in enumerate(cstags):
        if "I" in cs:
            cstags[i] = cs + cstags[i + 1][0]
            cstags[i + 1] = cstags[i + 1][1:]
    return cstags


def to_fixed_length(cs: str) -> str:
    if "D" in cs:
        cs = re.sub("[acgt]", "D,", cs[1:])
    elif "I" in cs:
        cs = f"{len(cs)-2}{cs[-1]},"
    elif "S" in cs:
        cs = "S,"
    elif "M" in cs:
        cs = cs.replace("M", "M,")
    return cs


def cstag_to_mids(cstag: str) -> str:
    cstags = format_cstag(cstag)
    cstags = append_next_mids_to_ins(cstags)
    cstags_fixlen = list(map(to_fixed_length, cstags))
    mids = "".join(cstags_fixlen)
    return mids


def padding(mids: str, pos: int, reflen: int) -> str:
    midslen = mids.count(",")
    left_pad = "=," * (int(pos) - 1)
    right_pad = "=," * (reflen - midslen - int(pos) + 1)
    mids_padding = "".join([left_pad, mids, right_pad])
    return mids_padding


def trim(mids_padding: str, reflen: int) -> str:
    mids_trim = ",".join(mids_padding.split(",")[0:reflen])
    return mids_trim


def mids_small_mutation(alignments_unique: list) -> list:
    read = alignments_unique[0]["alignments"]
    record = read.split("\t")
    samdict = dict(
        qname=record[0].replace(",", "_"),
        reflen=int(record[2]),
        pos=int(record[3]),
        cstag=[_ for _ in record if "cs:Z:" in _][0],
    )
    samdict["mids"] = cstag_to_mids(samdict["cstag"])
    mids_padding = padding(samdict["mids"], samdict["pos"], samdict["reflen"])
    mids_trim = trim(mids_padding, samdict["reflen"])
    output = ",".join([samdict["qname"], mids_trim]).rstrip(",")
    return output


def mids_large_mutation(alignments_duplicated: list) -> str:
    read_nums = len(alignments_duplicated)
    saminfo = list()
    for read in alignments_duplicated:
        record = read["alignments"].split("\t")
        samdict = dict(
            qname=record[0].replace(",", "_"),
            reflen=int(record[2]),
            pos=int(record[3]),
            cstag=[_ for _ in record if "cs:Z:" in _][0],
        )
        saminfo.append(samdict)
    saminfo = sorted(saminfo, key=lambda x: x["pos"])
    mids = [cstag_to_mids(s["cstag"]) for s in saminfo]
    _ = [saminfo[i].update({"mids": s}) for i, s in enumerate(mids)]
    # large deletion
    if read_nums == 2:
        left_len = saminfo[0]["mids"].count(",") - 1
        del_len = saminfo[1]["pos"] - saminfo[0]["pos"] - left_len
        del_seq = "D," * del_len
        mids_join = "".join([saminfo[0]["mids"], del_seq, saminfo[1]["mids"]])
    # large inversion
    elif read_nums == 3:
        midslow = saminfo[1]["mids"].lower()
        saminfo[1]["mids"] = midslow
        mids_join = "".join(
            [saminfo[0]["mids"], saminfo[1]["mids"], saminfo[2]["mids"]]
        )
    else:
        return ""
    mids_padding = padding(mids_join, saminfo[0]["pos"], saminfo[0]["reflen"])
    mids_trim = trim(mids_padding, saminfo[0]["reflen"])
    output = ",".join([saminfo[0]["qname"], mids_trim]).rstrip(",")
    return output


def to_mids(alignments: list) -> str:
    if len(alignments) == 1:
        output = mids_small_mutation(alignments)
    else:
        output = mids_large_mutation(alignments)
    return output


def sam_to_mids(sampath: str, threads: int) -> list:
    with open(sampath, "r") as f:
        sam = f.readlines()
    # SQ
    sqheaders = (s for s in sam if s.startswith("@SQ"))
    sqheaders = get_sqheaders(sqheaders)
    # Alignments
    alignments = (s for s in sam if not s.startswith("@"))
    alignments = (s for s in alignments if "cs:Z:" in s)
    alignments = [
        s.replace("\t" + s.split("\t")[2], "\t" + str(sqheaders[s.split("\t")[2]]))
        for s in alignments
    ]
    # Group by qname
    aligndict = [{"qname": s.split("\t")[0], "alignments": s} for s in alignments]
    aligndict = sorted(aligndict, key=lambda x: x["qname"])
    aligngroup = (list(group) for _, group in groupby(aligndict, lambda x: x["qname"]))
    with ProcessPoolExecutor(max_workers=threads) as executor:
        # MIDS conversion
        mids = list(executor.map(to_mids, aligngroup))
    return mids
