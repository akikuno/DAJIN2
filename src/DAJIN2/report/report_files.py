from __future__ import annotations
import textwrap


def to_fasta(header: str, cons_seq: str) -> str:
    header = ">" + header
    cons_seq_wrap = textwrap.wrap(cons_seq, 80)
    fasta = "\n".join([header, *cons_seq_wrap]) + "\n"
    return fasta
