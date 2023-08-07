from __future__ import annotations

import textwrap
import cstag
from pathlib import Path


def _to_fasta(header: str, cons_seq: str) -> str:
    header = ">" + header
    cons_seq_wrap = textwrap.wrap(cons_seq, 80)
    fasta = "\n".join([header, *cons_seq_wrap]) + "\n"
    return fasta


def to_fasta(TEMPDIR: Path | str, SAMPLE_NAME: str, cons_sequence: dict) -> None:
    for header, cons_seq in cons_sequence.items():
        path_output = Path(TEMPDIR, "report", "FASTA", SAMPLE_NAME, f"{SAMPLE_NAME}_{header}.fasta")
        path_output.write_text(_to_fasta(header, cons_seq))


def _to_html(header: str, cons_per: list[dict]) -> str:
    cons_cssplit = [max(cons, key=cons.get) for cons in cons_per]
    prev_cstag = cons_cssplit[0]
    cons_cstag = []
    cons_cstag.append(prev_cstag)
    for current_cstag in cons_cssplit[1:]:
        if "=" == prev_cstag[0] == current_cstag[0]:
            cons_cstag.append(current_cstag[1])
        elif "-" == prev_cstag[0] == current_cstag[0]:
            cons_cstag.append(current_cstag[1])
        elif "+" == current_cstag[0]:
            cs_ins = current_cstag.replace("+", "").split("|")
            ins = "+" + "".join(cs_ins[:-1]) + cs_ins[-1]
            cons_cstag.append(ins)
        else:
            cons_cstag.append(current_cstag)
        prev_cstag = current_cstag
    cons_cstag = "".join(cons_cstag)
    return cstag.to_html(cons_cstag, header)


def to_html(TEMPDIR: Path | str, SAMPLE_NAME: str, cons_percentage: dict) -> None:
    for header, cons_per in cons_percentage.items():
        path_output = Path(TEMPDIR, "report", "HTML", SAMPLE_NAME, f"{SAMPLE_NAME}_{header}.html")
        path_output.write_text(_to_html(header, cons_per))


def to_vcf(header: str, cons_per: list[dict]) -> str:
    pass
