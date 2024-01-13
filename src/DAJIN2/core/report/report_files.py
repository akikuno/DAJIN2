from __future__ import annotations

import textwrap
import cstag
from pathlib import Path
from DAJIN2.utils.cssplits_handler import convert_cssplits_to_cstag


def parse_fasta(file_path: Path | str) -> tuple[str, str]:
    """
    Parses a FASTA file and returns the header and concatenated sequence.

    :param file_path: Path to the FASTA file
    :return: A tuple with header (string) and sequence (string)
    """
    with open(file_path, "r") as f:
        lines = f.readlines()

    header = lines[0].strip().lstrip(">")
    sequence = "".join(line.strip() for line in lines[1:])

    return header, sequence


def _to_fasta(header: str, sequence: str) -> str:
    header = ">" + header
    sequence_wrapped = textwrap.wrap(sequence, 80)
    fasta = "\n".join([header, *sequence_wrapped]) + "\n"
    return fasta


def to_fasta(TEMPDIR: Path | str, SAMPLE_NAME: str, cons_sequence: dict) -> None:
    for header, sequence in cons_sequence.items():
        path_output = Path(TEMPDIR, "report", "FASTA", SAMPLE_NAME, f"{SAMPLE_NAME}_{header}.fasta")
        path_output.write_text(_to_fasta(f"{SAMPLE_NAME}_{header}", sequence))


def to_fasta_reference(TEMPDIR: Path | str, SAMPLE_NAME: str) -> None:
    for fasta in Path(TEMPDIR, SAMPLE_NAME, "fasta").glob("*.fasta"):
        header, sequence = parse_fasta(fasta)
        path_output = Path(TEMPDIR, "report", "FASTA", SAMPLE_NAME, f"{header}.fasta")
        path_output.write_text(_to_fasta(f"{SAMPLE_NAME}_{header}", sequence))


def _to_html(SAMPLE_NAME: str, header: str, cons_per: list[dict]) -> str:
    cons_cssplit = [max(cons, key=cons.get) for cons in cons_per]
    cons_cstag = convert_cssplits_to_cstag(cons_cssplit)
    return cstag.to_html(cons_cstag, f"{SAMPLE_NAME} {header.replace('_', ' ')}")


def to_html(TEMPDIR: Path | str, SAMPLE_NAME: str, cons_percentage: dict) -> None:
    for header, cons_per in cons_percentage.items():
        path_output = Path(TEMPDIR, "report", "HTML", SAMPLE_NAME, f"{SAMPLE_NAME}_{header}.html")
        path_output.write_text(_to_html(SAMPLE_NAME, header, cons_per))


def to_vcf(header: str, cons_per: list[dict]) -> str:
    pass
