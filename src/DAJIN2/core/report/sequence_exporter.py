from __future__ import annotations

import textwrap
from pathlib import Path

from DAJIN2.core.report.html_builder import to_html
from DAJIN2.utils import fileio


def convert_to_fasta(header: str, sequence: str) -> str:
    header = ">" + header
    sequence_wrapped = textwrap.wrap(sequence, 80)
    fasta = "\n".join([header, *sequence_wrapped]) + "\n"

    return fasta


def convert_to_html(
    TEMPDIR: Path, SAMPLE_NAME: str, FASTA_ALLELES: dict, allele: str, header: str, cons_midsv_tag: str
) -> str:
    path_midsv_sv = Path(TEMPDIR, SAMPLE_NAME, "midsv", f"consensus_{allele}_midsv.jsonl")
    is_sv_allele = False
    if path_midsv_sv.exists():
        is_sv_allele = True
        midsv_sv_allele = list(fileio.read_jsonl(path_midsv_sv))
    else:
        midsv_sv_allele = ["=" + base for base in list(FASTA_ALLELES[allele])]

    return to_html(
        midsv_sv_allele, cons_midsv_tag, allele, is_sv_allele, description=f"{SAMPLE_NAME} {header.replace('_', ' ')}"
    )


##################################################
# Export files
##################################################


def export_to_fasta(TEMPDIR: Path, SAMPLE_NAME: str, cons_sequence: dict) -> None:
    for key, sequence in cons_sequence.items():
        header = key.replace("|", "_")
        path_output = Path(TEMPDIR, "report", "FASTA", SAMPLE_NAME, f"{SAMPLE_NAME}_{header}.fasta")
        with open(path_output, "w", newline="\n", encoding="utf-8") as f:
            f.write(convert_to_fasta(f"{SAMPLE_NAME}_{header}", sequence))


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
        with open(path_output, "w", newline="\n", encoding="utf-8") as f:
            f.write(convert_to_fasta(f"{SAMPLE_NAME}_{header}", sequence))


def export_to_html(
    TEMPDIR: Path,
    SAMPLE_NAME: str,
    FASTA_ALLELES: dict[str, str],
    cons_midsv_tags: dict[str, str],
    map_name_allele: dict[str, str],
) -> None:
    for key, cons_midsv_tag in cons_midsv_tags.items():
        allele = map_name_allele[key]
        header = key.replace("|", "_")
        path_output = Path(TEMPDIR, "report", "HTML", SAMPLE_NAME, f"{SAMPLE_NAME}_{header}.html")
        with open(path_output, "w", newline="\n", encoding="utf-8") as f:
            f.write(convert_to_html(TEMPDIR, SAMPLE_NAME, FASTA_ALLELES, allele, header, cons_midsv_tag))
