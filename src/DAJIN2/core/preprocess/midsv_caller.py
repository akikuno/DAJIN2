from __future__ import annotations

import midsv

from pathlib import Path
from typing import Generator
from itertools import chain, groupby

from collections import Counter

from DAJIN2.utils import io, sam_handler, cssplits_handler


def has_inversion_in_splice(CIGAR: str) -> bool:
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
        if has_inversion_in_splice(cigar_splice):
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


def replace_internal_n_to_d(midsv_sample: Generator[list[dict]], sequence: str) -> Generator[list[dict]]:
    """
    Replace internal 'N's with 'D' in a given sequence.
    This function modifies the 'CSSPLIT' field in the input sample. It identifies
    the boundaries of consecutive 'N's and replaces any 'N' within these boundaries
    with the corresponding character from the provided sequence. 'N's at the boundaries
    remain unchanged.
    """
    for samp in midsv_sample:
        cssplits = samp["CSSPLIT"].split(",")

        left_idx_n, right_idx_n = cssplits_handler.find_n_boundaries(cssplits)

        for j, (cs, seq_char) in enumerate(zip(cssplits, sequence)):
            if left_idx_n < j < right_idx_n and cs == "N":
                cssplits[j] = f"-{seq_char}"

        samp["CSSPLIT"] = ",".join(cssplits)
        yield samp


def is_reverse_strand(flag: int) -> bool:
    """
    Determines if the read is mapped to the reverse strand.
    """
    # Check if the bit 4 (0-based) is set
    return bool(flag & 16)


def convert_flag_to_strand(midsv_sample: Generator[list[dict]]) -> Generator[list[dict]]:
    """Convert FLAG to STRAND (+ or -)"""
    for samp in midsv_sample:
        samp["STRAND"] = "-" if is_reverse_strand(int(samp["FLAG"])) else "+"
        del samp["FLAG"]
        yield samp


def filter_samples_by_n_proportion(midsv_sample: Generator[dict], threshold: int = 95) -> Generator[list[dict]]:
    """Filters out the samples from the input generator where the proportion of 'N' in the 'CSSPLIT' field is 95% or higher."""
    for samp in midsv_sample:
        cssplits = samp.get("CSSPLIT", "").split(",")
        count = Counter(cssplits)
        n_percentage = count["N"] / sum(count.values()) * 100
        if n_percentage < threshold:
            yield samp


###########################################################
# convert_consecutive_indels_to_match
###########################################################

"""
Due to alignment errors, there can be instances where a true match is mistakenly replaced with "insertion following a deletion".
For example, although it should be "=C,=T", it gets replaced by "-C,+C|=T". In such cases, a process is performed to revert it back to "=C,=T".
"""


def convert_consecutive_indels_to_match(cssplit: str) -> str:
    i = 0
    cssplit_reversed = cssplit.split(",")[::-1]
    while i < len(cssplit_reversed):
        current_cs = cssplit_reversed[i]

        if not current_cs.startswith("+"):
            i += 1
            continue

        insertions = [base.lstrip("+") for base in current_cs.split("|")[:-1]][::-1]

        # Extract deletions
        deletions = []

        for j in range(1, len(insertions) + 1):
            if i + j >= len(cssplit_reversed):
                break

            next_cs = cssplit_reversed[i + j]

            if not next_cs.startswith("-"):
                break

            deletions.append(next_cs.lstrip("-"))

        if insertions != deletions:
            i += 1
            continue

        # Format insertions
        cssplit_reversed[i] = current_cs.split("|")[-1]

        # Format deletions
        for k, insertion in enumerate(insertions, 1):
            cssplit_reversed[i + k] = cssplit_reversed[i + k].replace("-", "=")

        i += len(insertions) + 1

    return ",".join(cssplit_reversed[::-1])


def convert_consecutive_indels(midsv_sample: Generator) -> Generator[list[dict]]:
    for m in midsv_sample:
        m["CSSPLIT"] = convert_consecutive_indels_to_match(m["CSSPLIT"])
        yield m


def reallocate_insertions_within_deletion(midsv_sample: Generator) -> Generator[list[dict]]:
    for m in midsv_sample:
        m["CSSPLIT"] = cssplits_handler.reallocate_insertion_within_deletion(m["CSSPLIT"])
        yield m


###########################################################
# main
###########################################################


def generate_midsv(ARGS, is_control: bool = False, is_insertion: bool = False) -> None:
    name = ARGS.control_name if is_control else ARGS.sample_name

    for allele, sequence in ARGS.fasta_alleles.items():
        if Path(ARGS.tempdir, name, "midsv", f"{allele}.json").exists():
            continue

        if is_control and is_insertion:
            """
            Set the destination for midsv as `barcode01/midsv/insertion1_barcode02.json` when control is barcode01, sample is barcode02, and the allele is insertion1.
            """
            path_ont = Path(ARGS.tempdir, name, "sam", f"map-ont_{allele}_{ARGS.sample_name}.sam")
            path_splice = Path(ARGS.tempdir, name, "sam", f"splice_{allele}_{ARGS.sample_name}.sam")
            path_output_midsv = Path(ARGS.tempdir, name, "midsv", f"{allele}_{ARGS.sample_name}.json")
        else:
            """
            Set the destination for midsv as `barcode02/midsv/insertion1.json` when the sample is barcode02 and the allele is insertion1.
            """
            path_ont = Path(ARGS.tempdir, name, "sam", f"map-ont_{allele}.sam")
            path_splice = Path(ARGS.tempdir, name, "sam", f"splice_{allele}.sam")
            path_output_midsv = Path(ARGS.tempdir, name, "midsv", f"{allele}.json")

        sam_ont = sam_handler.remove_overlapped_reads(list(io.read_sam(path_ont)))
        sam_splice = sam_handler.remove_overlapped_reads(list(io.read_sam(path_splice)))
        qname_of_map_ont = extract_qname_of_map_ont(sam_ont, sam_splice)
        sam_of_map_ont = filter_sam_by_preset(sam_ont, qname_of_map_ont, preset="map-ont")
        sam_of_splice = filter_sam_by_preset(sam_splice, qname_of_map_ont, preset="splice")
        midsv_chaind = transform_to_midsv_format(chain(sam_of_map_ont, sam_of_splice))
        midsv_sample = replace_internal_n_to_d(midsv_chaind, sequence)
        midsv_sample = convert_flag_to_strand(midsv_sample)
        midsv_sample = filter_samples_by_n_proportion(midsv_sample)
        midsv_sample = convert_consecutive_indels(midsv_sample)
        midsv_sample = reallocate_insertions_within_deletion(midsv_sample)
        midsv.write_jsonl(midsv_sample, path_output_midsv)
