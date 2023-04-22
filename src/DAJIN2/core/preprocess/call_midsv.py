from pathlib import Path
from collections import defaultdict
import midsv

def _get_qname_of_single_read(path_sam: str | Path) -> set():
    count_qname = defaultdict(int)
    with open(path_sam) as f:
        for line in f:
            qname = line.split("\t")[0]
            count_qname[qname] += 1
    qname_of_single_read = set()
    for key, count in count_qname.items():
        if count == 1:
            qname_of_single_read.add(key)
    return qname_of_single_read

def _extract_sam(path_sam: str | Path, qname_of_single_read: set, single:bool=True) -> list:
    sam = midsv.read_sam(path_sam)
    sam_extracted = []
    for s in sam:
        if s[0].startswith("@"):
            sam_extracted.append(s)
        if single:
            if s[0] in qname_of_single_read:
                sam_extracted.append(s)
        else:
            if s[0] not in qname_of_single_read:
                sam_extracted.append(s)
    return sam_extracted

def _midsv_transform(sam: list[list[str]]) -> list[list[str]]:
    num_header = 0
    for s in sam:
        if s[0].startswith("@"):
            num_header += 1
        else:
            break
    if len(sam) == num_header:
        return []
    return midsv.transform(sam, midsv=False, cssplit=True, qscore=False)


def call_midsv(TEMPDIR, SAMPLE_NAME, allele) -> None:
    path_sam_ont = f"{TEMPDIR}/sam/{SAMPLE_NAME}_map-ont_{allele}.sam"
    path_sam_splice = f"{TEMPDIR}/sam/{SAMPLE_NAME}_splice_{allele}.sam"
    qname_of_single_read = _get_qname_of_single_read(path_sam_ont)
    sam_of_single_read = _extract_sam(path_sam_ont, qname_of_single_read)
    sam_of_multiple_reads = _extract_sam(path_sam_splice, qname_of_single_read, single=False)
    midsv_of_single_read = _midsv_transform(sam_of_single_read)
    midsv_of_multiple_reads = _midsv_transform(sam_of_multiple_reads)
    midsv_sample = midsv_of_single_read + midsv_of_multiple_reads
    midsv.write_jsonl(midsv_sample, f"{TEMPDIR}/midsv/{SAMPLE_NAME}_{allele}.jsonl")

# def call_midsv(TEMPDIR, path_sam):
#     sam = midsv.read_sam(path_sam)
#     midsv_jsonl = midsv.transform(sam, midsv=False, cssplit=True, qscore=True)
#     output_jsonl = Path(TEMPDIR, "midsv", f"{path_sam.stem}.jsonl")
#     midsv.write_jsonl(midsv_jsonl, output_jsonl)
