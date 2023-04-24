from pathlib import Path
from collections import defaultdict
import midsv
import re

def _get_qname_of_single_read(sam: list[list[str]]) -> set():
    count_qname = defaultdict(int)
    for alignment in sam:
        if alignment[0].startswith("@"):
            continue
        qname = alignment[0]
        count_qname[qname] += 1
    qname_of_single_read = set()
    for key, count in count_qname.items():
        if count == 1:
            qname_of_single_read.add(key)
    return qname_of_single_read


def _split_cigar(CIGAR:str) -> list[str]:
    cigar = re.split(r"([MIDNSH=X])", CIGAR)
    n = len(cigar)
    cigar_split = []
    for i, j in zip(range(0, n, 2), range(1, n, 2)):
        cigar_split.append(cigar[i] + cigar[j])
    return cigar_split

def _call_alignment_length(CIGAR: str) -> int:
    cigar_split = _split_cigar(CIGAR)
    alignment_length = 0
    for c in cigar_split:
        if re.search(r"[MDN=X]", c[-1]):
            alignment_length += int(c[:-1])
    return alignment_length

def _has_inversion_in_splice(CIGAR: str) -> bool:
    is_splice = False
    is_insertion = False
    for cigar in _split_cigar(CIGAR):
        if cigar.endswith("I"):
            is_insertion = True
            continue
        if is_insertion and cigar.endswith("N"):
            is_splice = True
            break
        else:
            is_insertion = False
    return is_splice


def _extract_qname_of_map_ont(sam_ont: list[list[str]], sam_splice: list[list[str]], qname_of_single_read: set) -> set():
    """Extract qname of reads from `map-ont` when:
        - single read
        - long alignment length
        - no inversion signal in `splice` alignment (insertion + deletion)
    """
    dict_sam_ont = {s[0]: s for s in sam_ont if s[0] in qname_of_single_read}
    dict_sam_splice = {s[0]: s for s in sam_splice if s[0] in qname_of_single_read}
    qname_of_map_ont = set()
    for qname in qname_of_single_read:
        if not qname in dict_sam_splice:
            continue
        alignment_ont = dict_sam_ont[qname]
        alignment_splice = dict_sam_splice[qname]
        if _has_inversion_in_splice(alignment_splice[5]):
            continue
        alignment_length_ont = _call_alignment_length(alignment_ont[5])
        alignment_length_splice = _call_alignment_length(alignment_splice[5])
        if alignment_length_ont >= alignment_length_splice:
            qname_of_map_ont.add(qname)
    return qname_of_map_ont


def _extract_sam(sam: list[list[str]], qname_of_map_ont: set, preset:str="map-ont") -> list[list[str]]:
    sam_extracted = []
    for alignment in sam:
        if alignment[0].startswith("@"):
            sam_extracted.append(alignment)
        if preset == "map-ont":
            if alignment[0] in qname_of_map_ont:
                sam_extracted.append(alignment)
        else:
            if alignment[0] not in qname_of_map_ont:
                sam_extracted.append(alignment)
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
    sam_ont = midsv.read_sam(f"{TEMPDIR}/sam/{SAMPLE_NAME}_map-ont_{allele}.sam")
    sam_splice = midsv.read_sam(f"{TEMPDIR}/sam/{SAMPLE_NAME}_splice_{allele}.sam")
    qname_of_single_read = _get_qname_of_single_read(sam_ont)
    qname_of_map_ont = _extract_qname_of_map_ont(sam_ont, sam_splice, qname_of_single_read)
    sam_of_map_ont = _extract_sam(sam_ont, qname_of_map_ont, preset="map-ont")
    sam_of_splice = _extract_sam(sam_splice, qname_of_map_ont, preset="splice")
    midsv_of_single_read = _midsv_transform(sam_of_map_ont)
    midsv_of_multiple_reads = _midsv_transform(sam_of_splice)
    midsv_sample = midsv_of_single_read + midsv_of_multiple_reads
    midsv.write_jsonl(midsv_sample, f"{TEMPDIR}/midsv/{SAMPLE_NAME}_{allele}.jsonl")

# def call_midsv(TEMPDIR, path_sam):
#     sam = midsv.read_sam(path_sam)
#     midsv_jsonl = midsv.transform(sam, midsv=False, cssplit=True, qscore=True)
#     output_jsonl = Path(TEMPDIR, "midsv", f"{path_sam.stem}.jsonl")
#     midsv.write_jsonl(midsv_jsonl, output_jsonl)
