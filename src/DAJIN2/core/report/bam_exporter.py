from __future__ import annotations

from collections import defaultdict
from itertools import groupby
from pathlib import Path

import pysam

from DAJIN2.utils import io, sam_handler


def recalculate_sam_coodinates_to_reference(sam: list[list[str]], GENOME_COODINATES: dict) -> list[str]:
    """Recalculate SAM genomic coordinates with the reference genome, not with the FASTA_ALLELE"""
    sam_headers = [s for s in sam if s[0].startswith("@")]
    sam_contents = [s for s in sam if not s[0].startswith("@")]
    for s in sam_headers:
        if s[0] != "@SQ":
            continue
        s[1] = f'SN:{GENOME_COODINATES["chrom"]}'
        s[2] = f'LN:{GENOME_COODINATES["chrom_size"]}'
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

    def convert_line(line: list[str]) -> list[str]:
        if not line[0].startswith("@") and line[3] == "0":
            line[3] = "1"
        return line

    return [convert_line(line) for line in sam_lines]


def group_by_name(sam_contents: list[str], clust_sample: list[dict]) -> dict[list]:
    """Group alignments in map-ont.sam by allele name (NAME)"""
    sam_contents.sort()
    clust_sample_sorted = sorted(clust_sample, key=lambda x: x["QNAME"])

    qnames: set[str] = {c["QNAME"] for c in clust_sample_sorted}

    sam_groups = defaultdict(list)
    idx_sam_contents = 0
    idx_clust_sample = 0
    while idx_sam_contents < len(sam_contents) and idx_clust_sample < len(clust_sample_sorted):
        alignments_sam = sam_contents[idx_sam_contents][:-1]  # Discard CS tags to reduce file size
        alignments_clsut_sample = clust_sample_sorted[idx_clust_sample]
        qname_sam = alignments_sam[0]

        if qname_sam not in qnames:
            idx_sam_contents += 1
            continue

        if qname_sam == alignments_clsut_sample["QNAME"]:
            key = alignments_clsut_sample["NAME"]
            sam_groups[key].append(alignments_sam)
            idx_sam_contents += 1
        else:
            idx_clust_sample += 1

    return dict(sam_groups)


###############################################################################
# igvjs
###############################################################################


def subset_qnames(RESULT_SAMPLE, readnum: int = 100) -> dict[set[str]]:
    qnames_by_name = defaultdict(set)
    for name, group in groupby(RESULT_SAMPLE, key=lambda x: x["NAME"]):
        group = list(group)
        qnames = [res["QNAME"] for res in group[:readnum]]
        qnames_by_name[name] = set(qnames)
    return dict(qnames_by_name)


def subset_reads(sam_content: list[str], qnames: set[str]) -> list[str]:
    return [sam for sam in sam_content if sam[0] in qnames]


###############################################################################
# Output
###############################################################################


def write_sam_to_bam(sam: list[list[str]], path_sam: str | Path, path_bam: str | Path, threads: int = 1) -> None:
    formatted_sam = "\n".join("\t".join(s) for s in sam)
    Path(path_sam).write_text(formatted_sam + "\n")
    pysam.sort("-@", f"{threads}", "-o", str(path_bam), str(path_sam))
    pysam.index("-@", f"{threads}", str(path_bam))


def update_sam(sam: list, GENOME_COODINATES: dict = None) -> list:
    if GENOME_COODINATES is None:
        GENOME_COODINATES = {}

    sam_records = sam.copy()
    sam_records = sam_handler.remove_microhomology(sam_records)
    if GENOME_COODINATES["genome"]:
        return recalculate_sam_coodinates_to_reference(sam_records, GENOME_COODINATES)
    else:
        return convert_pos_to_one_indexed(sam_records)


def export_to_bam(TEMPDIR, NAME, GENOME_COODINATES, THREADS, UUID, RESULT_SAMPLE=None, is_control=False) -> None:
    path_sam_input = Path(TEMPDIR, NAME, "sam", "control", "map-ont.sam")
    if not path_sam_input.exists():  # In the case of short-read.
        path_sam_input = Path(TEMPDIR, NAME, "sam", "control", "sr.sam")

    sam_records = list(io.read_sam(path_sam_input))

    # Update sam
    sam_updated = update_sam(sam_records, GENOME_COODINATES)

    # Output SAM and BAM
    path_sam_output = Path(TEMPDIR, "report", "BAM", f"temp_{UUID}_{NAME}_control.sam")
    path_bam_output = Path(TEMPDIR, "report", "BAM", NAME, f"{NAME}.bam")
    write_sam_to_bam(sam_updated, path_sam_output, path_bam_output, THREADS)

    # Prepare SAM headers and contents
    sam_headers = [s for s in sam_updated if s[0].startswith("@")]
    sam_contents = [s for s in sam_updated if not s[0].startswith("@")]
    if is_control:
        qnames_100reads: set[str] = set(list({s[0] for s in sam_contents[:10000]})[:100])  # subset 100 reads
        sam_subset = [s for s in sam_updated if s[0] in qnames_100reads]
        path_sam_output = Path(TEMPDIR, "report", "BAM", f"temp_{UUID}_{NAME}_control_cache.sam")
        path_bam_output = Path(TEMPDIR, "cache", ".igvjs", NAME, "control.bam")
        write_sam_to_bam(sam_headers + sam_subset, path_sam_output, path_bam_output, THREADS)
    else:
        sam_groups = group_by_name(sam_contents, RESULT_SAMPLE)
        qnames_by_name = subset_qnames(RESULT_SAMPLE)
        # Output SAM and BAM
        for name, sam_content in sam_groups.items():
            # BAM
            path_sam_output = Path(TEMPDIR, "report", "BAM", f"temp_{UUID}_{name}.sam")
            path_bam_output = Path(TEMPDIR, "report", "BAM", NAME, f"{NAME}_{name}.bam")
            write_sam_to_bam(sam_headers + sam_content, path_sam_output, path_bam_output, THREADS)
            # igvjs
            sam_subset = subset_reads(sam_content, qnames_by_name[name])
            path_sam_output = Path(TEMPDIR, "report", "BAM", f"temp_{UUID}_{name}_subset.sam")
            path_bam_output = Path(TEMPDIR, "report", ".igvjs", NAME, f"{name}.bam")
            write_sam_to_bam(sam_headers + sam_subset, path_sam_output, path_bam_output, THREADS)

    # Remove temporary files
    sam_temp = Path(TEMPDIR, "report", "BAM").glob(f"temp_{UUID}*.sam")
    [s.unlink() for s in sam_temp]
