from __future__ import annotations

import textwrap
from pathlib import Path

import cstag

from DAJIN2.core.report.insertion_reflector import reflect_ref_insertion_to_query
from DAJIN2.utils.cssplits_handler import convert_cssplits_to_cstag, reallocate_insertion_within_deletion


def convert_to_fasta(header: str, sequence: str) -> str:
    header = ">" + header
    sequence_wrapped = textwrap.wrap(sequence, 80)
    fasta = "\n".join([header, *sequence_wrapped]) + "\n"

    return fasta


def convert_to_html(TEMPDIR: Path, SAMPLE_NAME: str, header: str, cons_per: list[dict]) -> str:
    cons_cssplits = [max(cons, key=cons.get) for cons in cons_per]
    cons_cssplits_reallocated = reallocate_insertion_within_deletion(cons_cssplits, bin_size=500, percentage=50)
    cons_cstag = convert_cssplits_to_cstag(cons_cssplits_reallocated)

    allele = header.split("_")[1]
    if Path(TEMPDIR, SAMPLE_NAME, "cstag", f"{allele}.txt").exists():
        ref_cstag = Path(TEMPDIR, SAMPLE_NAME, "cstag", f"{allele}.txt").read_text()
        cons_cstag = reflect_ref_insertion_to_query(ref_cstag, cons_cstag)

    return cstag.to_html(cons_cstag, f"{SAMPLE_NAME} {header.replace('_', ' ')}")


##################################################
# Export files
##################################################


def export_to_fasta(TEMPDIR: Path, SAMPLE_NAME: str, cons_sequence: dict) -> None:
    for header, sequence in cons_sequence.items():
        path_output = Path(TEMPDIR, "report", "FASTA", SAMPLE_NAME, f"{SAMPLE_NAME}_{header}.fasta")
        path_output.write_text(convert_to_fasta(f"{SAMPLE_NAME}_{header}", sequence))


def parse_fasta(file_path: Path) -> tuple[str, str]:
    """Parses a FASTA file and returns the header and concatenated sequence."""
    with open(file_path) as f:
        lines = f.readlines()

    header = lines[0].strip().lstrip(">")
    sequence = "".join(line.strip() for line in lines[1:])

    return header, sequence


def export_reference_to_fasta(TEMPDIR: Path, SAMPLE_NAME: str) -> None:
    for fasta in Path(TEMPDIR, SAMPLE_NAME, "fasta").glob("*.fasta"):
        header, sequence = parse_fasta(fasta)
        path_output = Path(TEMPDIR, "report", "FASTA", SAMPLE_NAME, f"{header}.fasta")
        path_output.write_text(convert_to_fasta(f"{SAMPLE_NAME}_{header}", sequence))


def export_to_html(TEMPDIR: Path, SAMPLE_NAME: str, cons_percentage: dict[list]) -> None:
    for header, cons_per in cons_percentage.items():
        path_output = Path(TEMPDIR, "report", "HTML", SAMPLE_NAME, f"{SAMPLE_NAME}_{header}.html")
        path_output.write_text(convert_to_html(TEMPDIR, SAMPLE_NAME, header, cons_per))


# def to_vcf(TEMPDIR: Path, SAMPLE_NAME: str, GENOME_COODINATES: dict[str, str], cons_percentage: dict[list]) -> str:
#     pass
#     for header, cons_per in cons_percentage.items():
#         path_output = Path(TEMPDIR, "report", "HTML", SAMPLE_NAME, f"{SAMPLE_NAME}_{header}.vcf")
#         path_output.write_text(_to_html(TEMPDIR, SAMPLE_NAME, header, cons_per))
