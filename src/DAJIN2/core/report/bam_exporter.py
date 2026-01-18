from __future__ import annotations

import random
import uuid
from collections import defaultdict
from pathlib import Path

import pysam

from DAJIN2.core.preprocess.alignment.mapping import to_sam
from DAJIN2.utils import fileio, sam_handler


def recalculate_sam_coordinates_to_reference(sam: list[list[str]], genome_coordinates: dict) -> list[str]:
    """Recalculate SAM genomic coordinates with the reference genome, not with the FASTA_ALLELE"""
    sam_headers = [s for s in sam if s[0].startswith("@")]
    sam_contents = [s for s in sam if not s[0].startswith("@")]
    for s in sam_headers:
        if s[0] != "@SQ":
            continue
        s[1] = f"SN:{genome_coordinates['chrom']}"
        s[2] = f"LN:{genome_coordinates['chrom_size']}"

    for s in sam_contents:
        s[2] = genome_coordinates["chrom"]

    if genome_coordinates["strand"] == "-":
        sam_contents = sam_handler.revcomp_sam(sam_contents, genome_coordinates["end"])
    else:
        for s in sam_contents:
            s[3] = str(int(s[3]) + genome_coordinates["start"] - 1)

    return sam_headers + sam_contents


def convert_pos_to_one_indexed(sam_lines: list[list[str]]) -> list[list[str]]:
    """Convert SAM POS from 0-indexed to 1-indexed"""

    def convert_index(line: list[str]) -> list[str]:
        if not line[0].startswith("@") and line[3] == "0":
            line[3] = "1"
        return line

    return [convert_index(line) for line in sam_lines]


def group_by_allele_name(sam_contents: list[str], result_sample: list[dict]) -> dict[list]:
    """Group alignments in the input SAM by allele name (NAME)."""
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


def sample_unique_qnames(qnames: list[str], max_reads: int, rng: random.Random) -> set[str]:
    unique_qnames = list(dict.fromkeys(qnames))
    if len(unique_qnames) <= max_reads:
        return set(unique_qnames)
    return set(rng.sample(unique_qnames, max_reads))


def sample_sam_by_qname(sam_content: list[list[str]], max_reads: int, rng: random.Random) -> list[list[str]]:
    qnames = [sam[0] for sam in sam_content]
    sampled_qnames = sample_unique_qnames(qnames, max_reads, rng)
    return [sam for sam in sam_content if sam[0] in sampled_qnames]


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


def update_genome_coordinates(sam: list, genome_coordinates: dict | None) -> list:
    if genome_coordinates is None:
        genome_coordinates = {}

    sam_records = sam.copy()
    sam_records = sam_handler.remove_microhomology(sam_records)
    # Update @SQ headers if genome coordinates are available (from --genome or --genome-coordinate)
    if genome_coordinates.get("chrom") and genome_coordinates.get("chrom_size"):
        return recalculate_sam_coordinates_to_reference(sam_records, genome_coordinates)
    else:
        return convert_pos_to_one_indexed(sam_records)


def resolve_best_sam_path(tempdir: str | Path, name: str, allele: str = "control") -> Path:
    path_best_sam = Path(tempdir, name, "midsv", allele, f"{name}_best.sam")
    if path_best_sam.exists():
        return path_best_sam
    return Path(tempdir, name, "sam", allele, "map-ont.sam")


def export_sequence_error_to_bam(TEMPDIR, NAME, genome_coordinates, THREADS) -> None:
    path_fastq = Path(TEMPDIR, NAME, "fastq", f"{NAME}_sequence_error.fastq.gz")
    path_fasta = Path(TEMPDIR, NAME, "fasta", "control.fasta")
    preset = "map-ont"

    sam_records = [record.split("\t") for record in to_sam(path_fasta, path_fastq, preset=preset, threads=THREADS)]
    sam_updated = update_genome_coordinates(sam_records, genome_coordinates)

    # Output SAM and BAM
    path_bam_output = Path(TEMPDIR, "report", "BAM", NAME, f"{NAME}_sequence_error.bam")
    write_sam_to_bam(sam_updated, path_bam_output, THREADS)


def export_to_bam(TEMPDIR, NAME, genome_coordinates, THREADS, RESULT_SAMPLE=None, is_control=False) -> None:
    path_sam_input = resolve_best_sam_path(TEMPDIR, NAME, allele="control")

    sam_records = list(fileio.read_sam(path_sam_input))
    sam_updated = update_genome_coordinates(sam_records, genome_coordinates)

    rng = random.Random()
    max_reads_igvjs = 100

    # Output SAM and BAM
    path_bam_output = Path(TEMPDIR, "report", "BAM", NAME, f"{NAME}.bam")
    write_sam_to_bam(sam_updated, path_bam_output, THREADS)

    # Prepare SAM headers and contents
    sam_headers, sam_contents = split_headers_and_contents(sam_updated)

    if is_control:
        sam_subset = sample_sam_by_qname(sam_contents, max_reads_igvjs, rng)
        path_bam_output = Path(TEMPDIR, "cache", ".igvjs", NAME, "control.bam")
        path_bam_output.parent.mkdir(parents=True, exist_ok=True)
        write_sam_to_bam(sam_headers + sam_subset, path_bam_output, THREADS)
    else:
        sam_groups = group_by_allele_name(sam_contents, RESULT_SAMPLE)
        # Output SAM and BAM
        for name, sam_content in sam_groups.items():
            # BAM
            path_bam_output = Path(TEMPDIR, "report", "BAM", NAME, f"{NAME}_{name}.bam")
            write_sam_to_bam(sam_headers + sam_content, path_bam_output, THREADS)
            # igvjs
            sam_subset = sample_sam_by_qname(sam_content, max_reads_igvjs, rng)
            path_bam_output = Path(TEMPDIR, "report", ".igvjs", NAME, f"{name}.bam")
            write_sam_to_bam(sam_headers + sam_subset, path_bam_output, THREADS)

    # Export sequence error to BAM
    export_sequence_error_to_bam(TEMPDIR, NAME, genome_coordinates, THREADS)
