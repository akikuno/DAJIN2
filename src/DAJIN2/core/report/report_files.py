from __future__ import annotations
import cstag
import textwrap


def to_fasta(header: str, cons_seq: str) -> str:
    header = ">" + header
    cons_seq_wrap = textwrap.wrap(cons_seq, 80)
    fasta = "\n".join([header, *cons_seq_wrap]) + "\n"
    return fasta


def to_html(header: str, cons_per: list[dict]) -> str:
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


def to_vcf(header: str, cons_per: list[dict]) -> str:
    pass

