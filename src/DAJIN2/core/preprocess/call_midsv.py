from __future__ import annotations

import json
import re
from itertools import chain, groupby
from pathlib import Path
from typing import Generator

import midsv


def _load_sam(path_sam: str | Path) -> Generator[list[str]]:
    return midsv.read_sam(path_sam)


def _split_cigar(CIGAR: str) -> list[str]:
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


def _extract_qname_of_map_ont(sam_ont: Generator[list[str]], sam_splice: Generator[list[str]]) -> set():
    """Extract qname of reads from `map-ont` when:
    - no inversion signal in `splice` alignment (insertion + deletion)
    - single read
    - long alignment length
    """
    dict_alignments_splice = {s[0]: s for s in sam_splice if not s[0].startswith("@")}
    alignments_ont = [s for s in sam_ont if not s[0].startswith("@")]
    alignments_ont.sort(key=lambda x: x[0])
    qname_of_map_ont = set()
    for qname_ont, group in groupby(alignments_ont, key=lambda x: x[0]):
        alignment_ont = list(group)
        if qname_ont not in dict_alignments_splice:
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


def _extract_sam(sam: Generator[list[str]], qname_of_map_ont: set, preset: str = "map-ont") -> Generator[list[str]]:
    for alignment in sam:
        if alignment[0].startswith("@"):
            yield alignment
        if preset == "map-ont":
            if alignment[0] in qname_of_map_ont:
                yield alignment
        else:
            if alignment[0] not in qname_of_map_ont:
                yield alignment


def _midsv_transform(sam: Generator[list[str]]) -> Generator[list[str]]:
    for midsv_sample in midsv.transform(sam, midsv=False, cssplit=True, qscore=False):
        yield midsv_sample


def _replaceNtoD(midsv_sample: Generator[list[str]], sequence: str) -> Generator[dict[str, str]]:
    for samp in midsv_sample:
        qname = samp["QNAME"]
        cssplits = samp["CSSPLIT"].split(",")
        # extract right/left index of the end of sequential Ns
        left_idx_n = 0
        for cs in cssplits:
            if cs != "N":
                break
            left_idx_n += 1
        right_idx_n = 0
        for cs in cssplits[::-1]:
            if cs != "N":
                break
            right_idx_n += 1
        right_idx_n = len(cssplits) - right_idx_n - 1
        # replace sequential Ns within the sequence
        for j, (cs, seq) in enumerate(zip(cssplits, sequence)):
            if left_idx_n <= j <= right_idx_n and cs == "N":
                cssplits[j] = f"-{seq}"
        yield {"QNAME": qname, "CSSPLIT": ",".join(cssplits)}


def call_midsv(TEMPDIR: Path | str, FASTA_ALLELES: dict, SAMPLE_NAME: str) -> None:
    for allele, sequence in FASTA_ALLELES.items():
        path_ont = Path(TEMPDIR, "sam", f"{SAMPLE_NAME}_map-ont_{allele}.sam")
        path_splice = Path(TEMPDIR, "sam", f"{SAMPLE_NAME}_splice_{allele}.sam")
        qname_of_map_ont = _extract_qname_of_map_ont(_load_sam(path_ont), _load_sam(path_splice))
        sam_of_map_ont = _extract_sam(_load_sam(path_ont), qname_of_map_ont, preset="map-ont")
        sam_of_splice = _extract_sam(_load_sam(path_splice), qname_of_map_ont, preset="splice")
        sam_chained = chain(sam_of_map_ont, sam_of_splice)
        midsv_chaind = _midsv_transform(sam_chained)
        midsv_sample = _replaceNtoD(midsv_chaind, sequence)
        filepath = Path(TEMPDIR, "midsv", f"{SAMPLE_NAME}_{allele}.json")
        with open(filepath, "wt", encoding="utf-8") as f:
            for data in midsv_sample:
                f.write(json.dumps(data) + "\n")
