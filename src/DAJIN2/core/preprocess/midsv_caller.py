from __future__ import annotations

import midsv

from pathlib import Path
from typing import Generator
from itertools import chain, groupby

# from collections import Counter

from DAJIN2.utils import sam_handler
from DAJIN2.utils import cssplits_handler


def _has_inversion_in_splice(CIGAR: str) -> bool:
    previous_insertion = False
    for cigar in sam_handler.split_cigar(CIGAR):
        if cigar.endswith("I"):
            previous_insertion = True
            continue
        if previous_insertion and cigar.endswith("N"):
            return True
        else:
            previous_insertion = False
    return False


def extract_qname_of_map_ont(sam_ont: Generator[list[str]], sam_splice: Generator[list[str]]) -> set[str]:
    """Extract qname of reads from `map-ont` when:
    - no inversion signal in `splice` alignment (insertion + deletion)
    - single read
    - long alignment length
    """
    alignments_splice = {s[0]: s for s in sam_splice if not s[0].startswith("@")}
    alignments_ont = sorted([s for s in sam_ont if not s[0].startswith("@")], key=lambda x: x[0])
    qname_of_map_ont = set()
    for qname_ont, group in groupby(alignments_ont, key=lambda x: x[0]):
        group = list(group)

        if qname_ont not in alignments_splice:
            qname_of_map_ont.add(qname_ont)
            continue

        cigar_splice = alignments_splice[qname_ont][5]

        # If preset=splice and inversion is present, `midsv.transform`` will not work, so use preset=map-ont.
        if _has_inversion_in_splice(cigar_splice):
            qname_of_map_ont.add(qname_ont)
            continue

        # only accept single read or inversion reads
        # TODO: 逆位のリードは単純に「３つのリードがある」という条件でいいのか？
        if len(group) == 2 or len(group) >= 4:
            continue
        cigar_ont = group[0][5]

        alignment_length_ont = sam_handler.calculate_alignment_length(cigar_ont)
        alignment_length_splice = sam_handler.calculate_alignment_length(cigar_splice)

        if alignment_length_ont >= alignment_length_splice:
            qname_of_map_ont.add(qname_ont)
    return qname_of_map_ont


def filter_sam_by_preset(sam: Generator[list[str]], qname_of_map_ont: set, preset: str = "map-ont") -> Generator:
    for alignment in sam:
        if alignment[0].startswith("@"):
            yield alignment
        elif preset == "map-ont" and alignment[0] in qname_of_map_ont:
            yield alignment
        elif preset != "map-ont" and alignment[0] not in qname_of_map_ont:
            yield alignment


def transform_to_midsv_format(sam: Generator[list[str]]) -> Generator[list[dict]]:
    for midsv_sample in midsv.transform(sam, midsv=False, cssplit=True, qscore=False, keep=set(["FLAG"])):
        yield midsv_sample


def replace_n_to_d(midsv_sample: Generator[list[dict]], sequence: str) -> Generator[list[dict]]:
    """Replace N with D, but not contiguous from both ends."""
    for samp in midsv_sample:
        cssplits = samp["CSSPLIT"].split(",")

        left_idx_n, right_idx_n = cssplits_handler.find_n_boundaries(cssplits)

        # Replace Ns within the sequence boundaries with the corresponding sequence character
        for j, (cs, seq_char) in enumerate(zip(cssplits, sequence)):
            if left_idx_n < j < right_idx_n and cs == "N":
                cssplits[j] = f"-{seq_char}"

        samp["CSSPLIT"] = ",".join(cssplits)
        yield samp


def convert_flag_to_strand(midsv_sample: Generator[list[dict]]) -> Generator[list[dict]]:
    """Convert FLAG to STRAND (+ or -)"""
    REVERSE_STRAND_FLAG = 16 | 2064
    for samp in midsv_sample:
        samp["STRAND"] = "-" if samp["FLAG"] & REVERSE_STRAND_FLAG else "+"
        del samp["FLAG"]
        yield samp


# def filter_samples_by_n_proportion(midsv_sample: Generator[dict]) -> Generator[list[dict]]:
#     """Filters out the samples from the input generator where the proportion of 'N' in the 'CSSPLIT' field is 90% or higher."""
#     for m in midsv_sample:
#         cssplits = m.get("CSSPLIT", "").split(",")
#         n_count = Counter(cssplits)["N"]
#         total_count = len(cssplits)
#         if total_count == 0 or n_count / total_count > 0.9:
#             continue
#         else:
#             yield m


###########################################################
# main
###########################################################


def execute(TEMPDIR: Path | str, FASTA_ALLELES: dict, NAME: str) -> None:
    for allele, sequence in FASTA_ALLELES.items():
        path_output = Path(TEMPDIR, NAME, "midsv", f"{allele}.json")
        if path_output.exists():
            continue
        path_ont = Path(TEMPDIR, NAME, "sam", f"map-ont_{allele}.sam")
        path_splice = Path(TEMPDIR, NAME, "sam", f"splice_{allele}.sam")
        sam_ont = sam_handler.remove_overlapped_reads(list(midsv.read_sam(path_ont)))
        sam_splice = sam_handler.remove_overlapped_reads(list(midsv.read_sam(path_splice)))
        qname_of_map_ont = extract_qname_of_map_ont(sam_ont, sam_splice)
        sam_of_map_ont = filter_sam_by_preset(sam_ont, qname_of_map_ont, preset="map-ont")
        sam_of_splice = filter_sam_by_preset(sam_splice, qname_of_map_ont, preset="splice")
        midsv_chaind = transform_to_midsv_format(chain(sam_of_map_ont, sam_of_splice))
        midsv_sample = replace_n_to_d(midsv_chaind, sequence)
        midsv_sample = convert_flag_to_strand(midsv_sample)
        # midsv_sample = filter_samples_by_n_proportion(midsv_sample)
        midsv.write_jsonl(midsv_sample, path_output)
