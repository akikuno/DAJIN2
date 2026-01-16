from __future__ import annotations

import shutil
from collections import defaultdict
from collections.abc import Iterator
from pathlib import Path

import midsv

from DAJIN2.utils import cssplits_handler, fileio, sam_handler


def extract_preset_and_cigar_by_qname(path_sam_files: list[Path]) -> dict[dict[str, str]]:
    preset_cigar_by_qname = defaultdict(dict)
    # Extract preset and CIGAR
    for path in path_sam_files:
        preset = path.stem
        sam: list[list[str]] = fileio.read_sam(path)
        for record in sam:
            if record[0].startswith("@"):
                continue
            qname = record[0]
            cigar = record[5]
            if qname in preset_cigar_by_qname and preset in preset_cigar_by_qname[qname]:
                cigar += preset_cigar_by_qname[qname][preset]
            preset_cigar_by_qname[qname].update({preset: cigar})
    return dict(preset_cigar_by_qname)


def extract_best_preset(preset_cigar_by_qname: dict[str, dict[str, str]]) -> dict[str, str]:
    """Select the preset with the longest alignment length."""
    best_preset = defaultdict(str)
    for qname in preset_cigar_by_qname:
        preset_cigar = preset_cigar_by_qname[qname]
        if preset_cigar == {}:
            continue

        alignment_lengths = {
            preset: sam_handler.calculate_alignment_length(cigar) for preset, cigar in preset_cigar.items()
        }

        # If all alignment lengths are the same, prioritize map-ont
        if len(set(alignment_lengths.values())) == 1 and "map-ont" in alignment_lengths:
            best_preset[qname] = "map-ont"
            continue

        # Define a custom key function to prioritize map-ont
        def custom_key(key: str, alignment_lengths=alignment_lengths) -> tuple[int, bool]:
            return (alignment_lengths[key], key == "map-ont")

        max_key = max(alignment_lengths, key=custom_key)
        best_preset[qname] = max_key

    return dict(best_preset)


def extract_best_alignment_length_from_sam(path_sam_files: list[Path], best_preset: dict[str, str]) -> list[list[str]]:
    best_alignments: list[list[str]] = []
    flag_header = False
    for path in path_sam_files:
        preset = path.stem
        sam = fileio.read_sam(path)
        for record in sam:
            if record[0].startswith("@"):
                if not flag_header:
                    best_alignments.append(record)
            else:
                qname = record[0]
                if best_preset.get(qname) == preset:
                    best_alignments.append(record)
        flag_header = True
    return best_alignments


def write_sam(alignments: list[list[str]], path_output: Path) -> Path:
    path_output.parent.mkdir(parents=True, exist_ok=True)
    with open(path_output, "w") as f:
        for record in alignments:
            f.write("\t".join(record) + "\n")
    return path_output


def transform_to_midsv_format(path_sam: Path) -> list[dict]:
    return midsv.transform(path_sam=path_sam, qscore=False, keep={"FLAG"})


def replace_internal_n_to_d(midsv_sample: Iterator[list[dict]], sequence: str) -> Iterator[list[dict]]:
    """
    Replace internal 'N's with 'D' in a given sequence.
    This function modifies the 'MIDSV' field in the input sample. It identifies
    the boundaries of consecutive 'N's and replaces any 'N' within these boundaries
    with the corresponding character from the provided sequence. 'N's at the boundaries
    remain unchanged.
    """
    for samp in midsv_sample:
        midsv_tags = samp["MIDSV"].split(",")

        left_idx_n, right_idx_n = cssplits_handler.find_n_boundaries(midsv_tags)

        for j, (cs, seq_char) in enumerate(zip(midsv_tags, sequence)):
            if left_idx_n < j < right_idx_n and cssplits_handler.is_n_tag(cs):
                midsv_tags[j] = f"-{seq_char}"

        samp["MIDSV"] = ",".join(midsv_tags)
        yield samp


def is_reverse_strand(flag: int) -> bool:
    """
    Determines if the read is mapped to the reverse strand.
    """
    # Check if the bit 4 (0-based) is set
    return bool(flag & 16)


def convert_flag_to_strand(midsv_sample: Iterator[list[dict]]) -> Iterator[list[dict]]:
    """Convert FLAG to STRAND (+ or -)"""
    for samp in midsv_sample:
        samp["STRAND"] = "-" if is_reverse_strand(int(samp["FLAG"])) else "+"
        del samp["FLAG"]
        yield samp


def filter_samples_by_n_proportion(midsv_sample: Iterator[dict], threshold: int = 95) -> Iterator[list[dict]]:
    """Filters out the samples from the input Iterator where the proportion of 'N' in the 'MIDSV' field is 95% or higher."""
    for samp in midsv_sample:
        midsv_tags = samp.get("MIDSV", "").split(",")
        total = len(midsv_tags)
        if total == 0:
            continue
        n_count = sum(1 for tag in midsv_tags if cssplits_handler.is_n_tag(tag))
        n_percentage = n_count / total * 100
        if n_percentage < threshold:
            yield samp


###########################################################
# convert_consecutive_indels_to_match
###########################################################


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
        for k, _ in enumerate(insertions, 1):
            cssplit_reversed[i + k] = cssplit_reversed[i + k].replace("-", "=")

        i += len(insertions) + 1

    return ",".join(cssplit_reversed[::-1])


def convert_consecutive_indels(midsv_sample: Iterator) -> Iterator[list[dict]]:
    """
    Due to alignment errors, there can be instances where a true match is mistakenly replaced with "insertion following a deletion".
    For example, although it should be "=C,=T", it gets replaced by "-C,+C|=T". In such cases, a process is performed to revert it back to "=C,=T".
    """
    for m in midsv_sample:
        m["MIDSV"] = convert_consecutive_indels_to_match(m["MIDSV"])
        yield m


def append_best_preset(midsv_sample: list[dict], best_preset: dict[str, str]) -> Iterator[dict]:
    for m in midsv_sample:
        qname = m["QNAME"]
        if qname in best_preset:
            m["PRESET"] = best_preset[qname]
        yield m


###########################################################
# main
###########################################################


def generate_midsv(ARGS, is_control: bool = False, is_sv: bool = False) -> None:
    name = ARGS.control_name if is_control else ARGS.sample_name

    for allele, sequence in ARGS.fasta_alleles.items():
        # Skip if sam file does not exist
        if not Path(ARGS.tempdir, name, "sam", allele).exists():
            continue

        path_midsv_directory = Path(ARGS.tempdir, name, "midsv", allele)
        path_midsv_directory.mkdir(parents=True, exist_ok=True)
        path_output_jsonl = Path(path_midsv_directory, f"{name}_midsv.jsonl")
        # Skip if midsv jsonl already exists and is not empty
        if path_output_jsonl.exists() and path_output_jsonl.stat().st_size > 0:
            continue

        if is_control and is_sv:
            # Set the destination for midsv as `barcode01/midsv/insertion1/barcode02_midsv.jsonl`
            # when control is barcode01, sample is barcode02, and the allele is insertion1.
            path_sam_files = list(Path(ARGS.tempdir, name, "sam", allele).glob(f"{ARGS.sample_name}_*.sam"))
            path_midsv_output = Path(ARGS.tempdir, name, "midsv", allele, f"{ARGS.sample_name}_midsv.jsonl")
        else:
            # Set the destination for midsv as `barcode02/midsv/insertion1/barcode02_midsv.jsonl`
            # when the sample is barcode02 and the allele is insertion1.
            path_sam_files = list(Path(ARGS.tempdir, name, "sam", allele).glob("*.sam"))
            path_midsv_output = Path(ARGS.tempdir, name, "midsv", allele, f"{name}_midsv.jsonl")

        preset_cigar_by_qname = extract_preset_and_cigar_by_qname(path_sam_files)
        best_preset = extract_best_preset(preset_cigar_by_qname)
        sam_best_alignments = extract_best_alignment_length_from_sam(path_sam_files, best_preset)

        # Remove path_midsv_directory if there are no alignments
        if sum(1 for record in sam_best_alignments if not record[0].startswith("@")) == 0:
            shutil.rmtree(path_midsv_directory, ignore_errors=True)
            continue

        best_sam_name = f"{ARGS.sample_name}_best.sam" if is_control and is_sv else f"{name}_best.sam"
        path_best_sam = write_sam(sam_best_alignments, Path(path_midsv_directory, best_sam_name))
        midsv_chaind = transform_to_midsv_format(path_best_sam)
        midsv_sample = replace_internal_n_to_d(midsv_chaind, sequence)
        midsv_sample = convert_flag_to_strand(midsv_sample)
        midsv_sample = filter_samples_by_n_proportion(midsv_sample)
        midsv_sample = convert_consecutive_indels(midsv_sample)
        midsv_sample = append_best_preset(midsv_sample, best_preset)
        fileio.write_jsonl(midsv_sample, path_midsv_output)
