from __future__ import annotations

import midsv
import re
from itertools import groupby

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


def _extract_qname_of_map_ont(sam_ont: list[list[str]], sam_splice: list[list[str]]) -> set():
    """Extract qname of reads from `map-ont` when:
        - no inversion signal in `splice` alignment (insertion + deletion)
        - single read
        - long alignment length
    """
    alignments_ont = [s for s in sam_ont if not s[0].startswith("@")]
    alignments_ont.sort(key=lambda x: x[0])
    dict_alignments_splice = {s[0]: s for s in sam_splice if not s[0].startswith("@")}
    qname_of_map_ont = set()
    for qname_ont, group in groupby(alignments_ont, key=lambda x: x[0]):
        alignment_ont = list(group)
        if not qname_ont in dict_alignments_splice:
            qname_of_map_ont.add(qname_ont)
            continue
        alignment_splice = dict_alignments_splice[qname_ont]
        if _has_inversion_in_splice(alignment_splice[5]):
            qname_of_map_ont.add(qname_ont)
            continue
        if len(alignment_ont) != 1:
            continue
        alignment_ont = alignment_ont[0]
        alignment_length_ont = _call_alignment_length(alignment_ont[5])
        alignment_length_splice = _call_alignment_length(alignment_splice[5])
        if alignment_length_ont >= alignment_length_splice:
            qname_of_map_ont.add(qname_ont)
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
    qname_of_map_ont = _extract_qname_of_map_ont(sam_ont, sam_splice)
    sam_of_map_ont = _extract_sam(sam_ont, qname_of_map_ont, preset="map-ont")
    sam_of_splice = _extract_sam(sam_splice, qname_of_map_ont, preset="splice")
    midsv_of_single_read = _midsv_transform(sam_of_map_ont)
    midsv_of_multiple_reads = _midsv_transform(sam_of_splice)
    midsv_sample = midsv_of_single_read + midsv_of_multiple_reads
    midsv.write_jsonl(midsv_sample, f"{TEMPDIR}/midsv/{SAMPLE_NAME}_{allele}.jsonl")
