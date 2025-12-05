from __future__ import annotations

import re
import textwrap
from pathlib import Path

from DAJIN2.core.report.html_builder import to_html
from DAJIN2.utils import io


def extract_allele_from_header(header: str) -> str:
    """
    Extract the allele name from a header string.

    Header format: allele{ID}_{allele_name}_{suffix}_{percent}%
    Examples:
        - allele01_25003_Tombola_TMF2635-2636_intact_100% -> 25003_Tombola_TMF2635-2636
        - allele01_control_intact_100% -> control
        - allele02_deletion01_SV_75% -> deletion01

    Args:
        header (str): The header string to parse

    Returns:
        str: The extracted allele name
    """
    # Pattern to match: allele{digits}_{allele_name}_{suffix}_{percent}%
    # suffix can be: intact, indels, SV
    # percent can be integer or decimal
    pattern = r"^allele\d+_(.+)_(intact|indels|SV)_\d+(?:\.\d+)?%$"
    match = re.match(pattern, header)

    if match:
        return match.group(1)

    # Fallback to original method for unexpected formats
    return header.split("_")[1] if "_" in header else header


def convert_to_fasta(header: str, sequence: str) -> str:
    header = ">" + header
    sequence_wrapped = textwrap.wrap(sequence, 80)
    fasta = "\n".join([header, *sequence_wrapped]) + "\n"

    return fasta


def convert_to_html(
    TEMPDIR: Path,
    SAMPLE_NAME: str,
    FASTA_ALLELES: dict,
    header: str,
    cons_midsv_tag: list[str],
    sv_name_map: dict[str, str] | None = None,
) -> str:
    allele = extract_allele_from_header(header)
    sv_name_map = sv_name_map or {}
    display_to_internal = {display: internal for internal, display in sv_name_map.items()}
    allele_internal = display_to_internal.get(allele, allele)

    path_midsv_sv = Path(TEMPDIR, SAMPLE_NAME, "midsv", f"consensus_{allele_internal}.jsonl")
    is_sv_allele = False
    if path_midsv_sv.exists():
        is_sv_allele = True
        midsv_sv_allele = list(io.read_jsonl(path_midsv_sv))
    else:
        allele_key = allele_internal if allele_internal in FASTA_ALLELES else allele
        midsv_sv_allele = ["=" + base for base in list(FASTA_ALLELES[allele_key])]

    return to_html(
        midsv_sv_allele, cons_midsv_tag, allele, is_sv_allele, description=f"{SAMPLE_NAME} {header.replace('_', ' ')}"
    )


##################################################
# Export files
##################################################


def export_to_fasta(TEMPDIR: Path, SAMPLE_NAME: str, cons_sequence: dict) -> None:
    for header, sequence in cons_sequence.items():
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


def export_reference_to_fasta(TEMPDIR: Path, SAMPLE_NAME: str, sv_name_map: dict[str, str] | None = None) -> None:
    sv_name_map = sv_name_map or {}
    for fasta in Path(TEMPDIR, SAMPLE_NAME, "fasta").glob("*.fasta"):
        header, sequence = parse_fasta(fasta)
        display_header = sv_name_map.get(header, header)
        path_output = Path(TEMPDIR, "report", "FASTA", SAMPLE_NAME, f"{display_header}.fasta")
        path_output.parent.mkdir(parents=True, exist_ok=True)
        with open(path_output, "w", newline="\n", encoding="utf-8") as f:
            f.write(convert_to_fasta(f"{SAMPLE_NAME}_{display_header}", sequence))


def export_to_html(
    TEMPDIR: Path, SAMPLE_NAME: str, FASTA_ALLELES: dict, cons_midsv_tags: dict[list], sv_name_map: dict[str, str] | None = None
) -> None:
    for header, cons_midsv_tag in cons_midsv_tags.items():
        path_output = Path(TEMPDIR, "report", "HTML", SAMPLE_NAME, f"{SAMPLE_NAME}_{header}.html")
        with open(path_output, "w", newline="\n", encoding="utf-8") as f:
            f.write(convert_to_html(TEMPDIR, SAMPLE_NAME, FASTA_ALLELES, header, cons_midsv_tag, sv_name_map))


# TODO: Implement to_vcf

# def to_vcf(TEMPDIR: Path, SAMPLE_NAME: str, GENOME_COODINATES: dict[str, str], cons_percentage: dict[list]) -> str:
#     pass
#     for header, cons_per in cons_percentage.items():
#         path_output = Path(TEMPDIR, "report", "HTML", SAMPLE_NAME, f"{SAMPLE_NAME}_{header}.vcf")
#         path_output.write_text(_to_html(TEMPDIR, SAMPLE_NAME, header, cons_per))
