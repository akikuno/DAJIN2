# from ctypes import alignment
import re
from itertools import groupby
from concurrent.futures import ProcessPoolExecutor

###############################################################################
# MIDS Conversion
###############################################################################


def extract_name_length(sam: list) -> dict:
    """
    Extract SN (Reference sequence name) and LN (Reference sequence length) information at SQ header from SAM file
    """
    sqheaders = (s for s in sam if s.startswith("@SQ"))
    SNLN = {}
    for sqheader in sqheaders:
        sn_ln = [sq for sq in sqheader.split("\t") if re.search(("SN:|LN:"), sq)]
        sn = sn_ln[0].replace("SN:", "")
        ln = sn_ln[1].replace("LN:", "")
        SNLN.update({sn: ln})
    return SNLN


def slide_insertion(CSTAGS: list) -> list:
    """
    Input:  ['MMMM', 'S', 'M', 'Dg', 'M', 'It', 'MMMM']
    Output: ['MMMM', 'S', 'M', 'Dg', 'M', 'ItM', 'MMM']
    """
    for i, cs in enumerate(CSTAGS):
        if "I" in cs:
            CSTAGS[i] = cs + CSTAGS[i + 1][0]
            CSTAGS[i + 1] = CSTAGS[i + 1][1:]
    return CSTAGS


def to_fixed_length(CSTAG: str) -> str:
    if "D" in CSTAG:
        CSTAG = re.sub("[acgt]", "D,", CSTAG[1:])
    elif "I" in CSTAG:
        CSTAG = f"{len(CSTAG)-2}{CSTAG[-1]},"
    elif "S" in CSTAG:
        CSTAG = "S,"
    elif "M" in CSTAG:
        CSTAG = CSTAG.replace("M", "M,")
    return CSTAG


def cstag_to_mids(CSTAG: str) -> str:
    """
    Input:  "cs:Z:=ACGT*ag=C-g=T+t=ACGT"
    Output: "M,M,M,M,S,M,D,M,1M,M,M,M"
    """
    # CSTAG = "cs:Z:=ACGT*ag=C-g=T+t=ACGT"
    CSTAG = CSTAG.replace("cs:Z:", "")
    CSTAG = CSTAG.replace("-", "=D")
    CSTAG = CSTAG.replace("+", "=I")
    CSTAG = re.sub("[ACGT]", "M", CSTAG)
    CSTAG = re.sub("\\*[acgt][acgt]", "=S", CSTAG)
    CSTAGS = CSTAG.split("=")[1:]
    CSTAGS = slide_insertion(CSTAGS)
    CSTAGS_fixlen = [to_fixed_length(cs) for cs in CSTAGS]
    return "".join(CSTAGS_fixlen).rstrip(",")


def padding(mids: str, pos: int, reflen: int) -> str:
    midslen = mids.count(",") + 1
    left_pad = "=," * (int(pos) - 1)
    right_pad = ",=" * (reflen - midslen - int(pos) + 1)
    return "".join([left_pad, mids, right_pad])


def trim(mids_padding: str, reflen: int) -> str:
    """
    Trim bases that are longer than the reference sequence
    """
    return ",".join(mids_padding.split(",")[0:reflen])


def mids_small_mutation(alignment: list) -> list:
    # alignment = aligngroupby[2]
    record = alignment[0]["alignment"].split("\t")
    idx = [i for i, a in enumerate(record) if "cs:Z" in a][0]
    samdict = dict(
        qname=record[0].replace(",", "_"),
        reflen=int(record[-1]),
        pos=int(record[3]),
        qual=record[10],
        cstag=record[idx],
    )
    samdict["mids"] = cstag_to_mids(samdict["cstag"])
    mids_padding = padding(samdict["mids"], samdict["pos"], samdict["reflen"])
    mids_trim = trim(mids_padding, samdict["reflen"])
    return ",".join([samdict["qname"], mids_trim])


def mids_large_deletion(alignments: list) -> str:
    saminfo = list()
    for read in alignments:
        record = read["alignment"].split("\t")
        idx = [i for i, a in enumerate(record) if "cs:Z" in a][0]
        samdict = dict(
            qname=record[0].replace(",", "_"),
            reflen=int(record[-1]),
            pos=int(record[3]),
            qual=record[10],
            cstag=record[idx],
        )
        saminfo.append(samdict)
    saminfo = sorted(saminfo, key=lambda x: x["pos"])
    mids = [cstag_to_mids(s["cstag"]) for s in saminfo]
    _ = [saminfo[i].update({"mids": s}) for i, s in enumerate(mids)]
    left_len = saminfo[0]["mids"].count(",") - 1
    del_len = saminfo[1]["pos"] - saminfo[0]["pos"] - left_len
    del_seq = "D," * del_len
    mids_join = "".join([saminfo[0]["mids"], del_seq, saminfo[1]["mids"]])
    mids_padding = padding(mids_join, saminfo[0]["pos"], saminfo[0]["reflen"])
    mids_trim = trim(mids_padding, saminfo[0]["reflen"])
    return ",".join([saminfo[0]["qname"], mids_trim])


def mids_large_inversion(alignments: list) -> str:
    saminfo = list()
    for read in alignments:
        record = read["alignment"].split("\t")
        idx = [i for i, a in enumerate(record) if "cs:Z" in a][0]
        samdict = dict(
            qname=record[0].replace(",", "_"),
            reflen=int(record[-1]),
            pos=int(record[3]),
            qual=record[10],
            cstag=record[idx],
        )
        saminfo.append(samdict)
    saminfo = sorted(saminfo, key=lambda x: x["pos"])
    mids = [cstag_to_mids(s["cstag"]) for s in saminfo]
    _ = [saminfo[i].update({"mids": s}) for i, s in enumerate(mids)]
    midslow = saminfo[1]["mids"].lower()
    saminfo[1]["mids"] = midslow
    mids_join = "".join([saminfo[0]["mids"], saminfo[1]["mids"], saminfo[2]["mids"]])
    mids_padding = padding(mids_join, saminfo[0]["pos"], saminfo[0]["reflen"])
    mids_trim = trim(mids_padding, saminfo[0]["reflen"])
    return ",".join([saminfo[0]["qname"], mids_trim])


def to_mids(aligngroupby: list) -> str:
    if len(aligngroupby) == 1:
        output = mids_small_mutation(aligngroupby)
    elif len(aligngroupby) == 2:
        output = mids_large_deletion(aligngroupby)
    elif len(aligngroupby) == 3:
        output = mids_large_inversion(aligngroupby)
    else:
        output = ""
    return output


def sam_to_mids(sampath: str, threads: int) -> list:
    with open(sampath, "r") as f:
        sam = f.read().splitlines()
    # SQ
    sqheaders = extract_name_length(sam)
    # Alignments
    alignments = []
    for alignment in sam:
        if "cs:Z:" not in alignment:
            continue
        # Avoid reads with too long SoftClip
        CIGAR = alignment.split("\t")[5]
        SOFTCLIPS = sum(int(S[:-1]) for S in re.split("([0-9]+S)", CIGAR) if "S" in S)
        SEQLEN = len(alignment.split("\t")[11])
        if SOFTCLIPS > SEQLEN / 2:
            continue
        RNAME = alignment.split("\t")[2]
        alignments.append("\t".join([alignment, sqheaders[RNAME]]))
    # Group by QNAME
    aligndict = [{"QNAME": a.split("\t")[0], "alignment": a} for a in alignments]
    aligndict = sorted(aligndict, key=lambda x: x["QNAME"])
    aligngroupby = [list(group) for _, group in groupby(aligndict, lambda x: x["QNAME"])]
    with ProcessPoolExecutor(max_workers=threads) as executor:
        # MIDS conversion
        mids = list(executor.map(to_mids, aligngroupby))
    return mids


###############################################################################
# Extract full-length leads only
###############################################################################


def extract_full_length_reads(mids: str) -> str:
    """
    Extract full-length leads only
    """
    mids = mids.split(",")
    left_cut = mids[1:51].count("=")
    right_cut = mids[-50:].count("=")
    if len(mids) > 100 and left_cut < 50 and right_cut < 50:
        return ",".join(mids)


###############################################################################
# MEMO
###############################################################################
# to_mids(aligngroupby[2])

# import collections
# c = collections.Counter()
# len(aligngroupby)
# len3 = []
# for a in aligngroupby:
#     c[len(a)] += 1
#     if len(a) == 3:
#         len3.append(a[0]["QNAME"])

# print(*len3, sep="\n")
# c

# sampath = ".tmpDAJIN/sam/barcode31_control.sam"

# with ProcessPoolExecutor(max_workers=threads) as executor:
#     # MIDS conversion
#     mids = list(executor.map(to_mids, aligngroupby))

# tmp_mids = mids[0:5]
