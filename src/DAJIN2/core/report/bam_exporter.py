from __future__ import annotations

import uuid
from collections import defaultdict
from itertools import groupby
from pathlib import Path

import pysam

from DAJIN2.core.preprocess.mapping import to_sam
from DAJIN2.utils import io, sam_handler


def recalculate_sam_coodinates_to_reference(sam: list[list[str]], GENOME_COODINATES: dict) -> list[str]:
    """Recalculate SAM genomic coordinates with the reference genome, not with the FASTA_ALLELE"""
    sam_headers = [s for s in sam if s[0].startswith("@")]
    sam_contents = [s for s in sam if not s[0].startswith("@")]
    for s in sam_headers:
        if s[0] != "@SQ":
            continue
        s[1] = f"SN:{GENOME_COODINATES['chrom']}"
        s[2] = f"LN:{GENOME_COODINATES['chrom_size']}"

    for s in sam_contents:
        s[2] = GENOME_COODINATES["chrom"]

    if GENOME_COODINATES["strand"] == "-":
        sam_contents = sam_handler.revcomp_sam(sam_contents, GENOME_COODINATES["end"])
    else:
        for s in sam_contents:
            s[3] = str(int(s[3]) + GENOME_COODINATES["start"] - 1)

    return sam_headers + sam_contents


def convert_pos_to_one_indexed(sam_lines: list[list[str]]) -> list[list[str]]:
    """Convert SAM POS from 0-indexed to 1-indexed"""

    def convert_index(line: list[str]) -> list[str]:
        if not line[0].startswith("@") and line[3] == "0":
            line[3] = "1"
        return line

    return [convert_index(line) for line in sam_lines]


def group_by_allele_name(sam_contents: list[str], result_sample: list[dict]) -> dict[list]:
    """Group alignments in map-ont.sam by allele name (NAME)"""
    sam_contents.sort()
    result_sample_sorted = sorted(result_sample, key=lambda x: x["QNAME"])

    qnames: set[str] = {c["QNAME"] for c in result_sample_sorted}

    sam_groups = defaultdict(list)
    idx_sam_contents = 0
    idx_result_sample = 0
    while idx_sam_contents < len(sam_contents) and idx_result_sample < len(result_sample_sorted):
        alignments_sam = sam_contents[idx_sam_contents][:-1]  # Discard CS tags to reduce file size
        alignments_clsut_sample = result_sample_sorted[idx_result_sample]
        qname_sam = alignments_sam[0]

        if qname_sam not in qnames:
            idx_sam_contents += 1
            continue

        if qname_sam == alignments_clsut_sample["QNAME"]:
            key = alignments_clsut_sample["NAME"]
            sam_groups[key].append(alignments_sam)
            idx_sam_contents += 1
        else:
            idx_result_sample += 1

    return dict(sam_groups)


###############################################################################
# igvjs
###############################################################################


def subset_qnames(result_sample: list[dict], readnum: int = 100) -> dict[set[str]]:
    qnames_by_name = defaultdict(set)
    for name, group in groupby(result_sample, key=lambda x: x["NAME"]):
        group = list(group)
        qnames = [res["QNAME"] for res in group[:readnum]]
        qnames_by_name[name] = set(qnames)
    return dict(qnames_by_name)


def subset_reads(sam_content: list[str], qnames: set[str]) -> list[str]:
    return [sam for sam in sam_content if sam[0] in qnames]


###############################################################################
# Output
###############################################################################


def split_headers_and_contents(sam: list) -> tuple:
    """Split SAM headers and contents."""
    sam_headers = [s for s in sam if s[0].startswith("@")]
    sam_contents = [s for s in sam if not s[0].startswith("@")]
    return sam_headers, sam_contents


def write_sam_to_bam(sam: list[list[str]], path_bam: str | Path, threads: int = 1) -> None:
    # Format SAM
    formatted_sam = "\n".join("\t".join(s) for s in sam) + "\n"

    # Write SAM to a temporary file
    path_sam = Path(path_bam.parent, f"temp_{uuid.uuid4()}.sam")
    with open(path_sam, "w", newline="\n", encoding="utf-8") as f:
        f.write(formatted_sam)

    # Convert SAM to BAM
    pysam.sort("-@", f"{threads}", "-o", str(path_bam), str(path_sam))
    pysam.index("-@", f"{threads}", str(path_bam))
    # Remove temporary file
    path_sam.unlink()


def update_genome_coodinates(sam: list, GENOME_COODINATES: dict = None) -> list:
    if GENOME_COODINATES is None:
        GENOME_COODINATES = {}

    sam_records = sam.copy()
    sam_records = sam_handler.remove_microhomology(sam_records)
    if GENOME_COODINATES["genome"]:
        return recalculate_sam_coodinates_to_reference(sam_records, GENOME_COODINATES)
    else:
        return convert_pos_to_one_indexed(sam_records)


def export_sequence_error_to_bam(TEMPDIR, NAME, GENOME_COODINATES, THREADS) -> None:
    path_fastq = Path(TEMPDIR, NAME, "fastq", f"{NAME}_sequence_error.fastq.gz")
    path_fasta = Path(TEMPDIR, NAME, "fasta", "control.fasta")
    preset = "map-ont"

    sam_records = [record.split("\t") for record in to_sam(path_fasta, path_fastq, preset=preset, threads=THREADS)]
    sam_updated = update_genome_coodinates(sam_records, GENOME_COODINATES)

    # Output SAM and BAM
    path_bam_output = Path(TEMPDIR, "report", "BAM", NAME, f"{NAME}_sequence_error.bam")
    write_sam_to_bam(sam_updated, path_bam_output, THREADS)


def export_to_bam(TEMPDIR, NAME, GENOME_COODINATES, THREADS, RESULT_SAMPLE=None, is_control=False) -> None:
    path_sam_input = Path(TEMPDIR, NAME, "sam", "control", "map-ont.sam")

    sam_records = list(io.read_sam(path_sam_input))
    sam_updated = update_genome_coodinates(sam_records, GENOME_COODINATES)

    # Output SAM and BAM
    path_bam_output = Path(TEMPDIR, "report", "BAM", NAME, f"{NAME}.bam")
    write_sam_to_bam(sam_updated, path_bam_output, THREADS)

    # Prepare SAM headers and contents
    sam_headers, sam_contents = split_headers_and_contents(sam_updated)

    if is_control:
        # subset 100 reads with uqniue qnames for igv.js
        qnames_100reads: set[str] = set(list({s[0] for s in sam_contents[:10000]})[:100])
        sam_subset = [s for s in sam_updated if s[0] in qnames_100reads]
        path_bam_output = Path(TEMPDIR, "cache", ".igvjs", NAME, "control.bam")
        write_sam_to_bam(sam_headers + sam_subset, path_bam_output, THREADS)
    else:
        sam_groups = group_by_allele_name(sam_contents, RESULT_SAMPLE)
        qnames_by_name = subset_qnames(RESULT_SAMPLE)
        # Output SAM and BAM
        for name, sam_content in sam_groups.items():
            # BAM
            path_bam_output = Path(TEMPDIR, "report", "BAM", NAME, f"{NAME}_{name}.bam")
            write_sam_to_bam(sam_headers + sam_content, path_bam_output, THREADS)
            # igvjs
            sam_subset = subset_reads(sam_content, qnames_by_name[name])
            path_bam_output = Path(TEMPDIR, "report", ".igvjs", NAME, f"{name}.bam")
            write_sam_to_bam(sam_headers + sam_subset, path_bam_output, THREADS)

    # Export sequence error to BAM
    export_sequence_error_to_bam(TEMPDIR, NAME, GENOME_COODINATES, THREADS)
