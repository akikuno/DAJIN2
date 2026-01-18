from __future__ import annotations

from pathlib import Path

from DAJIN2.utils.midsv_handler import revcomp_midsvs

LARGE_SV_THRESHOLD = 50


def _format_info(info: dict[str, str | int]) -> str:
    if not info:
        return "."
    priority = ["SVTYPE", "TYPE", "SVLEN", "SEQ", "QNAME"]
    seen = set()
    info_items: list[tuple[str, str | int]] = []
    for key in priority:
        if key in info:
            info_items.append((key, info[key]))
            seen.add(key)
    for key in sorted(k for k in info.keys() if k not in seen):
        info_items.append((key, info[key]))
    return ";".join(f"{key}={value}" for key, value in info_items)


def _reference_base(token: str) -> str:
    if token.startswith("+"):
        anchor = next((p for p in reversed(token.split("|")) if not p.startswith("+")), "")
        token = anchor or "=N"
    if len(token) < 2:
        return "N"
    ref_base = token[1]
    return ref_base.upper()


def _is_inversion_token(token: str) -> bool:
    return any(char.islower() for char in token if char.isalpha())


def _parse_insertion_token(token: str) -> tuple[str, str, str, str]:
    parts = token.split("|")
    inserted = "".join(part[1:].upper() for part in parts if part.startswith("+"))
    anchor = next((part for part in reversed(parts) if not part.startswith("+")), "")
    anchor_op = anchor[0] if anchor else "="
    anchor_body = anchor[1:].upper() if len(anchor) > 1 else "N"
    anchor_ref = anchor_body[0] if anchor_body else "N"
    anchor_alt = anchor_ref
    if anchor_op == "*" and len(anchor_body) >= 2:
        anchor_alt = anchor_body[1:]
    elif anchor_op == "-":
        anchor_alt = ""
    return inserted, anchor_ref, anchor_alt, anchor_op


def _is_unknown_token(token: str) -> bool:
    upper = token.upper()
    return upper in {"=N", "N"}


def _midsv_to_vcf_records(
    midsv_tags: list[str], chrom: str, start_offset: int, allele_id: str, large_sv_threshold: int
) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    pos = start_offset + 1
    inv_start = None
    inv_bases: list[str] = []
    n_start = None
    n_len = 0

    def flush_inversion() -> None:
        nonlocal inv_start, inv_bases
        if inv_start is None:
            return
        ref_seq = "".join(inv_bases).upper() or "N"
        info = {"SVTYPE": "INV", "SVLEN": len(inv_bases), "SEQ": ref_seq, "QNAME": allele_id}
        records.append(
            {"CHROM": chrom, "POS": inv_start, "ID": allele_id, "REF": ref_seq[0], "ALT": "<INV>", "INFO": info}
        )
        inv_start = None
        inv_bases = []

    def flush_unknown_run() -> None:
        nonlocal n_start, n_len
        if n_start is None:
            return
        ref_seq = "N" * n_len
        info = {"TYPE": "DEL", "SVLEN": -n_len, "SEQ": ref_seq, "QNAME": allele_id}
        records.append(
            {"CHROM": chrom, "POS": n_start, "ID": allele_id, "REF": ref_seq[0], "ALT": "<DEL>", "INFO": info}
        )
        n_start = None
        n_len = 0

    idx = 0
    while idx < len(midsv_tags):
        token = midsv_tags[idx]

        if _is_inversion_token(token):
            flush_unknown_run()
            if inv_start is None:
                inv_start = pos
            inv_bases.append(_reference_base(token))
            pos += 1
            idx += 1
            continue
        flush_inversion()

        if _is_unknown_token(token):
            if n_start is None:
                n_start = pos
                n_len = 0
            n_len += 1
            pos += 1
            idx += 1
            continue
        flush_unknown_run()

        if token.startswith("="):
            pos += 1
            idx += 1
            continue

        if token.startswith("*"):
            ref_base = _reference_base(token)
            alt_base = token[2:].upper() if len(token) >= 3 else ref_base
            info = {"TYPE": "SUB", "QNAME": allele_id}
            records.append(
                {"CHROM": chrom, "POS": pos, "ID": allele_id, "REF": ref_base, "ALT": alt_base, "INFO": info}
            )
            pos += 1
            idx += 1
            continue

        if token.startswith("-"):
            start_pos = pos
            deleted_seq = token[1:].upper()
            consumed = 1
            for look_ahead in midsv_tags[idx + 1 :]:
                if not look_ahead.startswith("-") or _is_inversion_token(look_ahead):
                    break
                deleted_seq += look_ahead[1:].upper()
                consumed += 1
            ref_base = deleted_seq[0] if deleted_seq else "N"
            info = {"TYPE": "DEL", "SVLEN": -len(deleted_seq), "SEQ": deleted_seq, "QNAME": allele_id}
            records.append(
                {"CHROM": chrom, "POS": start_pos, "ID": allele_id, "REF": ref_base, "ALT": "<DEL>", "INFO": info}
            )
            pos += consumed
            idx += consumed
            continue

        if token.startswith("+"):
            inserted_seq, anchor_ref, anchor_alt, anchor_op = _parse_insertion_token(token)
            info = {"TYPE": "INS", "SVLEN": len(inserted_seq), "QNAME": allele_id}
            if inserted_seq:
                info["SEQ"] = inserted_seq
            alt = "<INS>" if len(inserted_seq) > large_sv_threshold else anchor_ref + inserted_seq
            records.append(
                {"CHROM": chrom, "POS": pos, "ID": allele_id, "REF": anchor_ref, "ALT": alt, "INFO": info}
            )
            if anchor_op == "*" and anchor_alt and anchor_alt != anchor_ref:
                sub_info = {"TYPE": "SUB", "QNAME": allele_id}
                records.append(
                    {"CHROM": chrom, "POS": pos, "ID": allele_id, "REF": anchor_ref, "ALT": anchor_alt, "INFO": sub_info}
                )
            pos += 1
            idx += 1
            continue

        pos += 1
        idx += 1

    flush_inversion()
    flush_unknown_run()

    return records


def export_to_vcf(
    tempdir: Path, sample_name: str, genome_coordinates: dict | None, cons_midsv_tags: dict[str, list[str]]
) -> None:
    if genome_coordinates is None:
        genome_coordinates = {}

    chrom = str(genome_coordinates.get("chrom") or "control")
    start_offset = int(genome_coordinates.get("start", 0))

    def write_vcf(path_output: Path, records: list[dict[str, object]]) -> None:
        for order, record in enumerate(records):
            record["_order"] = order

        records.sort(key=lambda r: (r["CHROM"], r["POS"], r["_order"]))

        path_output.parent.mkdir(parents=True, exist_ok=True)
        with open(path_output, "w", newline="\n", encoding="utf-8") as f:
            f.write("##fileformat=VCFv4.3\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for record in records:
                info = record.get("INFO", {})
                if isinstance(info, dict) and "END" not in info:
                    svtype = info.get("SVTYPE") or info.get("TYPE")
                    if svtype == "DEL":
                        svlen = info.get("SVLEN")
                        try:
                            svlen_int = int(svlen)
                        except (TypeError, ValueError):
                            svlen_int = None
                        if svlen_int is not None:
                            info["END"] = int(record["POS"]) + abs(svlen_int) - 1
                info_str = _format_info(info)
                f.write(
                    f"{record['CHROM']}\t{record['POS']}\t{record.get('ID', '.')}\t{record['REF']}\t"
                    f"{record['ALT']}\t.\tPASS\t{info_str}\n"
                )

    for key, cons_midsv_tag in cons_midsv_tags.items():
        header = key.replace("|", "_")
        midsv_tags = cons_midsv_tag
        if genome_coordinates.get("strand") == "-":
            midsv_tags = revcomp_midsvs(midsv_tags)
        records = _midsv_to_vcf_records(midsv_tags, chrom, start_offset, header, LARGE_SV_THRESHOLD)
        path_output = Path(tempdir, "report", "VCF", sample_name, f"{sample_name}_{header}.vcf")
        write_vcf(path_output, records)
