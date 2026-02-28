from __future__ import annotations

from pathlib import Path
from tempfile import NamedTemporaryFile

import midsv

from DAJIN2.core.preprocess.alignment import mapping
from DAJIN2.utils import fileio
from DAJIN2.utils.allele_handler import to_allele_key
from DAJIN2.utils.dna_handler import revcomp
from DAJIN2.utils.midsv_handler import convert_midsvs_to_sequence
from DAJIN2.utils.report_name_handler import build_report_filename

LARGE_SV_THRESHOLD = 50
VALID_DNA_BASES = frozenset("ACGTN")
INFO_FIELD_PRIORITY = ("SVTYPE", "TYPE", "SVLEN", "SEQ", "QNAME")
MAPPY_PRESET = "map-ont"


def _build_vcf_record(
    chromosome_name: str,
    position: int,
    record_id: str,
    reference_bases: str,
    alternate_bases: str,
    info_fields: dict[str, str | int],
) -> dict[str, object]:
    """Create one internal VCF record object."""
    return {
        "CHROM": chromosome_name,
        "POS": position,
        "ID": record_id,
        "REF": reference_bases,
        "ALT": alternate_bases,
        "INFO": info_fields,
    }


def _format_info(info_fields: dict[str, str | int]) -> str:
    """Serialize INFO fields in a stable order for VCF output."""
    if not info_fields:
        return "."

    ordered_pairs: list[tuple[str, str | int]] = []
    emitted_keys: set[str] = set()

    for key in INFO_FIELD_PRIORITY:
        if key in info_fields:
            ordered_pairs.append((key, info_fields[key]))
            emitted_keys.add(key)

    for key in sorted(key for key in info_fields if key not in emitted_keys):
        ordered_pairs.append((key, info_fields[key]))

    serialized_pairs: list[str] = []
    for key, value in ordered_pairs:
        if isinstance(value, str):
            value = value.replace("%", "%25")
        serialized_pairs.append(f"{key}={value}")

    return ";".join(serialized_pairs)


def _reference_base(midsv_token: str) -> str:
    """Return the reference base implied by a MIDSV token."""
    token_for_reference = midsv_token
    if midsv_token.startswith("+"):
        anchor_token = next(
            (part for part in reversed(midsv_token.split("|")) if not part.startswith("+")),
            "",
        )
        token_for_reference = anchor_token or "=N"

    if len(token_for_reference) < 2:
        return "N"

    return token_for_reference[1].upper()


def _normalize_variant_sequence(sequence: str) -> str:
    """Strip non-DNA symbols from a sequence before VCF conversion."""
    return "".join(base for base in str(sequence).upper() if base in VALID_DNA_BASES)


def _revcomp_variant_sequence(sequence: str) -> str:
    """Reverse complement a sequence after removing MIDSV operators."""
    sanitized_sequence = _normalize_variant_sequence(sequence)
    return revcomp(sanitized_sequence) if sanitized_sequence else ""


def _is_inversion_token(midsv_token: str) -> bool:
    """Detect whether a MIDSV token contains inversion bases."""
    return any(character.islower() for character in midsv_token if character.isalpha())


def _parse_insertion_token(midsv_token: str) -> tuple[str, str, str, str]:
    """Extract inserted bases and anchor context from an insertion token."""
    token_parts = midsv_token.split("|")
    inserted_sequence = "".join(part[1:].upper() for part in token_parts if part.startswith("+"))
    anchor_token = next((part for part in reversed(token_parts) if not part.startswith("+")), "")
    anchor_operation = anchor_token[0] if anchor_token else "="
    anchor_body = anchor_token[1:].upper() if len(anchor_token) > 1 else "N"
    anchor_reference_base = anchor_body[0] if anchor_body else "N"
    anchor_alternate_bases = anchor_reference_base

    if anchor_operation == "*" and len(anchor_body) >= 2:
        anchor_alternate_bases = anchor_body[1:]
    elif anchor_operation == "-":
        anchor_alternate_bases = ""

    return anchor_operation, anchor_reference_base, anchor_alternate_bases, inserted_sequence


def _is_unknown_token(midsv_token: str) -> bool:
    """Detect placeholder tokens that represent unknown bases."""
    return midsv_token.upper() in {"=N", "N"}


def _build_inversion_record(
    chromosome_name: str,
    start_position: int,
    inversion_bases: list[str],
    allele_record_id: str,
) -> dict[str, object] | None:
    """Build a symbolic inversion record from buffered inversion bases."""
    if start_position is None:
        return None

    reference_sequence = "".join(inversion_bases).upper() or "N"
    info_fields: dict[str, str | int] = {
        "SVTYPE": "INV",
        "SVLEN": len(inversion_bases),
        "SEQ": reference_sequence,
        "QNAME": allele_record_id,
    }
    return _build_vcf_record(
        chromosome_name,
        start_position,
        allele_record_id,
        reference_sequence[0],
        "<INV>",
        info_fields,
    )


def _build_unknown_run_record(
    chromosome_name: str,
    start_position: int,
    run_length: int,
    allele_record_id: str,
) -> dict[str, object] | None:
    """Represent an unknown-base run as a symbolic deletion-like record."""
    if start_position is None:
        return None

    reference_sequence = "N" * run_length
    info_fields: dict[str, str | int] = {
        "SVTYPE": "DEL",
        "TYPE": "DEL",
        "SVLEN": -run_length,
        "SEQ": reference_sequence,
        "QNAME": allele_record_id,
    }
    return _build_vcf_record(
        chromosome_name,
        start_position,
        allele_record_id,
        reference_sequence[0],
        "<DEL>",
        info_fields,
    )


def _collect_deleted_sequence(midsv_tokens: list[str], start_index: int) -> tuple[str, int]:
    """Merge consecutive deletion tokens into one deleted sequence."""
    deleted_bases = midsv_tokens[start_index][1:].upper()
    consumed_token_count = 1

    for lookahead_token in midsv_tokens[start_index + 1 :]:
        if not lookahead_token.startswith("-") or _is_inversion_token(lookahead_token):
            break
        deleted_bases += lookahead_token[1:].upper()
        consumed_token_count += 1

    return deleted_bases, consumed_token_count


def _build_insertion_alt(reference_base: str, inserted_sequence: str, large_sv_threshold: int) -> str:
    """Choose symbolic or explicit ALT notation for an insertion."""
    if len(inserted_sequence) > large_sv_threshold:
        return "<INS>"
    return reference_base + inserted_sequence


def _midsv_to_vcf_records(
    midsv_tags: list[str],
    chrom: str,
    start_offset: int,
    allele_id: str,
    large_sv_threshold: int,
) -> list[dict[str, object]]:
    """Convert MIDSV tokens into local-coordinate VCF records."""
    records: list[dict[str, object]] = []
    local_position = start_offset + 1
    left_anchor_base: str | None = None
    inversion_start_position: int | None = None
    inversion_bases: list[str] = []
    unknown_run_start_position: int | None = None
    unknown_run_length = 0

    def flush_inversion_record() -> None:
        """Emit the buffered inversion event, if present."""
        nonlocal inversion_start_position, inversion_bases
        inversion_record = _build_inversion_record(chrom, inversion_start_position, inversion_bases, allele_id)
        if inversion_record is not None:
            records.append(inversion_record)
        inversion_start_position = None
        inversion_bases = []

    def flush_unknown_run_record() -> None:
        """Emit the buffered unknown-base run, if present."""
        nonlocal unknown_run_start_position, unknown_run_length
        unknown_run_record = _build_unknown_run_record(
            chrom,
            unknown_run_start_position,
            unknown_run_length,
            allele_id,
        )
        if unknown_run_record is not None:
            records.append(unknown_run_record)
        unknown_run_start_position = None
        unknown_run_length = 0

    # Walk through MIDSV tokens from left to right while tracking the current
    # local reference position and any event that spans multiple tokens.
    token_index = 0
    while token_index < len(midsv_tags):
        current_token = midsv_tags[token_index]

        # Inversions may cover consecutive lowercase tokens, so keep buffering
        # them until a non-inversion token appears.
        if _is_inversion_token(current_token):
            flush_unknown_run_record()
            if inversion_start_position is None:
                inversion_start_position = local_position
            current_reference_base = _reference_base(current_token)
            inversion_bases.append(current_reference_base)
            left_anchor_base = current_reference_base
            local_position += 1
            token_index += 1
            continue
        flush_inversion_record()

        # Unknown-base runs are also merged into one symbolic event.
        if _is_unknown_token(current_token):
            if unknown_run_start_position is None:
                unknown_run_start_position = local_position
            unknown_run_length += 1
            left_anchor_base = _reference_base(current_token)
            local_position += 1
            token_index += 1
            continue
        flush_unknown_run_record()

        # Pure matches only advance the reference cursor and refresh the last
        # known anchor base for downstream insertions.
        if current_token.startswith("="):
            left_anchor_base = _reference_base(current_token)
            local_position += 1
            token_index += 1
            continue

        # Substitutions become one single-base VCF record at the current
        # reference position.
        if current_token.startswith("*"):
            reference_base = _reference_base(current_token)
            alternate_base = current_token[2:].upper() if len(current_token) >= 3 else reference_base
            substitution_info: dict[str, str | int] = {"TYPE": "SUB", "QNAME": allele_id}
            records.append(
                _build_vcf_record(
                    chrom,
                    local_position,
                    allele_id,
                    reference_base,
                    alternate_base,
                    substitution_info,
                )
            )
            left_anchor_base = reference_base
            local_position += 1
            token_index += 1
            continue

        # Consecutive deletions are collapsed into one VCF deletion record with
        # a computed END coordinate.
        if current_token.startswith("-"):
            deleted_sequence, consumed_token_count = _collect_deleted_sequence(midsv_tags, token_index)
            deletion_length = len(deleted_sequence)
            deletion_end_position = local_position + deletion_length - 1
            deletion_info: dict[str, str | int] = {
                "SVTYPE": "DEL",
                "TYPE": "DEL",
                "SVLEN": -deletion_length,
                "SEQ": deleted_sequence,
                "QNAME": allele_id,
                "END": deletion_end_position,
            }
            records.append(
                _build_vcf_record(
                    chrom,
                    local_position,
                    allele_id,
                    deleted_sequence[0] if deleted_sequence else "N",
                    "<DEL>",
                    deletion_info,
                )
            )
            if deleted_sequence:
                left_anchor_base = deleted_sequence[-1]
            local_position += consumed_token_count
            token_index += consumed_token_count
            continue

        # Insertions use the previous matched base as the left anchor when
        # available. If the anchor token itself encodes a substitution, emit an
        # additional substitution record at the same locus.
        if current_token.startswith("+"):
            (
                anchor_operation,
                anchor_reference_base,
                anchor_alternate_bases,
                inserted_sequence,
            ) = _parse_insertion_token(current_token)
            has_left_anchor = left_anchor_base is not None and local_position > start_offset + 1
            insertion_position = local_position - 1 if has_left_anchor else local_position
            insertion_reference_base = left_anchor_base if has_left_anchor else anchor_reference_base
            insertion_info: dict[str, str | int] = {"TYPE": "INS", "SVLEN": len(inserted_sequence), "QNAME": allele_id}
            if inserted_sequence:
                insertion_info["SEQ"] = inserted_sequence
            records.append(
                _build_vcf_record(
                    chrom,
                    insertion_position,
                    allele_id,
                    insertion_reference_base,
                    _build_insertion_alt(insertion_reference_base, inserted_sequence, large_sv_threshold),
                    insertion_info,
                )
            )
            if anchor_operation == "*" and anchor_alternate_bases and anchor_alternate_bases != anchor_reference_base:
                substitution_info: dict[str, str | int] = {"TYPE": "SUB", "QNAME": allele_id}
                records.append(
                    _build_vcf_record(
                        chrom,
                        local_position,
                        allele_id,
                        anchor_reference_base,
                        anchor_alternate_bases,
                        substitution_info,
                    )
                )
            left_anchor_base = anchor_reference_base
            local_position += 1
            token_index += 1
            continue

        # Any other token shape is treated conservatively as one consumed
        # reference base so the parser stays synchronized.
        left_anchor_base = _reference_base(current_token)
        local_position += 1
        token_index += 1

    # Flush buffered multi-token events that reach the end of the sequence.
    flush_inversion_record()
    flush_unknown_run_record()
    return records


def _should_align_to_control(allele_name: str) -> bool:
    """Return whether an allele should be re-aligned against control."""
    return allele_name.lower() != "control"


def _extract_allele_name(consensus_key: str) -> str:
    """Extract the allele name portion from a consensus header."""
    key_parts = consensus_key.split("|")
    if len(key_parts) >= 2:
        return key_parts[1]
    return consensus_key


def _load_control_sequence(tempdir: Path, sample_name: str) -> str | None:
    """Load the control FASTA sequence used for VCF alignment."""
    path_control_fasta = Path(tempdir, sample_name, "fasta", f"{to_allele_key('control')}.fasta")
    if not path_control_fasta.exists():
        return None

    for fasta_record in fileio.read_fasta(path_control_fasta):
        return fasta_record["sequence"]
    return None


def _write_temp_fasta(sequence: str, record_name: str) -> Path:
    """Write a temporary FASTA file and return its path."""
    with NamedTemporaryFile("w", delete=False, suffix=".fasta") as temporary_fasta:
        temporary_fasta.write(f">{record_name}\n")
        for start_index in range(0, len(sequence), 80):
            temporary_fasta.write(sequence[start_index : start_index + 80] + "\n")
        return Path(temporary_fasta.name)


def _filter_primary_sam_lines(sam_lines: list[str]) -> list[str]:
    """Remove supplementary alignments before MIDSV conversion."""
    filtered_lines: list[str] = []
    for sam_line in sam_lines:
        normalized_line = sam_line if sam_line.endswith("\n") else f"{sam_line}\n"
        if normalized_line.startswith("@"):
            filtered_lines.append(normalized_line)
            continue

        sam_fields = normalized_line.rstrip("\n").split("\t")
        if len(sam_fields) > 1:
            try:
                alignment_flag = int(sam_fields[1])
            except ValueError:
                alignment_flag = 0
            if alignment_flag & 2048:
                continue

        filtered_lines.append(normalized_line)

    return filtered_lines


def _align_sequence_to_control(control_sequence: str, query_sequence: str) -> list[str]:
    """Align a consensus sequence to control and return MIDSV tokens."""
    if control_sequence.upper() == query_sequence.upper():
        return [f"={base}" for base in control_sequence.upper()]

    path_reference_fasta: Path | None = None
    path_query_fasta: Path | None = None
    path_alignment_sam: Path | None = None

    try:
        path_reference_fasta = _write_temp_fasta(control_sequence, "control")
        path_query_fasta = _write_temp_fasta(query_sequence, "query")
        sam_lines = list(mapping.to_sam(path_reference_fasta, path_query_fasta, preset=MAPPY_PRESET))
        if not sam_lines:
            return []

        primary_sam_lines = _filter_primary_sam_lines(sam_lines)
        with NamedTemporaryFile("w", delete=False, suffix=".sam") as temporary_sam:
            path_alignment_sam = Path(temporary_sam.name)
            temporary_sam.writelines(primary_sam_lines)

        midsv_alignments = midsv.transform(path_sam=path_alignment_sam, qscore=False)
        if not midsv_alignments:
            return []

        return midsv_alignments[0]["MIDSV"].split(",")
    finally:
        for temporary_path in (path_reference_fasta, path_query_fasta, path_alignment_sam):
            if temporary_path is not None:
                temporary_path.unlink(missing_ok=True)


def _convert_record_on_negative_strand(
    record: dict[str, object],
    info_fields: dict[str, str | int],
    region_end: int,
) -> dict[str, object]:
    """Convert one local-coordinate record to genomic coordinates on the minus strand."""
    variant_kind = str(info_fields.get("TYPE") or info_fields.get("SVTYPE") or "")
    local_position = int(record["POS"])
    local_end_position = int(info_fields.get("END", local_position))
    converted_record = dict(record)
    sequence_in_info = str(info_fields.get("SEQ", ""))

    if variant_kind == "INS":
        genomic_position = region_end - local_position + 1
        reference_base = _revcomp_variant_sequence(str(record["REF"])) or "N"
        inserted_sequence = _revcomp_variant_sequence(sequence_in_info)
        alternate_bases = "<INS>" if str(record["ALT"]) == "<INS>" else reference_base + inserted_sequence
        if inserted_sequence:
            info_fields["SEQ"] = inserted_sequence
        converted_record.update({"POS": genomic_position, "REF": reference_base, "ALT": alternate_bases})
        return converted_record

    if variant_kind == "DEL":
        genomic_position = region_end - local_end_position + 1
        genomic_end_position = region_end - local_position + 1
        deleted_sequence = _revcomp_variant_sequence(sequence_in_info)
        reference_base = (
            deleted_sequence[0] if deleted_sequence else _revcomp_variant_sequence(str(record["REF"])) or "N"
        )
        info_fields["SEQ"] = deleted_sequence
        info_fields["END"] = genomic_end_position
        converted_record.update({"POS": genomic_position, "REF": reference_base, "ALT": "<DEL>"})
        return converted_record

    if variant_kind == "SUB":
        genomic_position = region_end - local_position + 1
        converted_record.update(
            {
                "POS": genomic_position,
                "REF": _revcomp_variant_sequence(str(record["REF"])) or "N",
                "ALT": _revcomp_variant_sequence(str(record["ALT"])) or "N",
            }
        )
        return converted_record

    if variant_kind == "INV":
        inverted_sequence = _revcomp_variant_sequence(sequence_in_info)
        genomic_position = region_end - local_end_position + 1
        info_fields["SEQ"] = inverted_sequence
        info_fields["END"] = region_end - local_position + 1
        converted_record.update(
            {
                "POS": genomic_position,
                "REF": inverted_sequence[0]
                if inverted_sequence
                else _revcomp_variant_sequence(str(record["REF"])) or "N",
                "ALT": "<INV>",
            }
        )
        return converted_record

    genomic_position = region_end - local_position + 1
    converted_record["POS"] = genomic_position
    if isinstance(record.get("REF"), str):
        converted_record["REF"] = _revcomp_variant_sequence(str(record["REF"])) or "N"
    if isinstance(record.get("ALT"), str) and not str(record["ALT"]).startswith("<"):
        converted_record["ALT"] = _revcomp_variant_sequence(str(record["ALT"])) or "N"
    return converted_record


def _convert_records_to_genomic_coordinates(
    records: list[dict[str, object]],
    genome_coordinates: dict[str, object] | None,
) -> list[dict[str, object]]:
    """Translate local VCF records into genomic coordinates."""
    if not records or genome_coordinates is None:
        return records

    chromosome_name = str(genome_coordinates.get("chrom") or "control")
    region_start = int(genome_coordinates.get("start", 1))
    region_end = int(genome_coordinates.get("end", region_start))
    strand = str(genome_coordinates.get("strand") or "+")

    genomic_records: list[dict[str, object]] = []
    for record in records:
        info_fields = dict(record.get("INFO", {})) if isinstance(record.get("INFO", {}), dict) else {}
        genomic_record = dict(record)
        genomic_record["CHROM"] = chromosome_name

        if strand != "-":
            genomic_record["POS"] = region_start + int(record["POS"]) - 1
            if "END" in info_fields:
                info_fields["END"] = region_start + int(info_fields["END"]) - 1
        else:
            genomic_record = _convert_record_on_negative_strand(genomic_record, info_fields, region_end)
            genomic_record["CHROM"] = chromosome_name

        genomic_record["INFO"] = info_fields
        genomic_records.append(genomic_record)

    genomic_records.sort(key=lambda record: (str(record["CHROM"]), int(record["POS"]), str(record["ID"])))
    return genomic_records


def _ensure_deletion_end(record: dict[str, object]) -> None:
    """Populate END for deletion records when it can be derived from SVLEN."""
    info_fields = record.get("INFO", {})
    if not isinstance(info_fields, dict) or "END" in info_fields:
        return

    variant_kind = info_fields.get("SVTYPE") or info_fields.get("TYPE")
    if variant_kind != "DEL":
        return

    sv_length = info_fields.get("SVLEN")
    try:
        deletion_length = abs(int(sv_length))
    except (TypeError, ValueError):
        return

    info_fields["END"] = int(record["POS"]) + deletion_length - 1


def _write_vcf(path_output: Path, records: list[dict[str, object]]) -> None:
    """Write internal VCF records to a VCF file."""
    sortable_records: list[dict[str, object]] = []
    for output_order, record in enumerate(records):
        record_with_order = dict(record)
        record_with_order["_order"] = output_order
        sortable_records.append(record_with_order)

    sortable_records.sort(key=lambda record: (str(record["CHROM"]), int(record["POS"]), int(record["_order"])))

    path_output.parent.mkdir(parents=True, exist_ok=True)
    with open(path_output, "w", newline="\n", encoding="utf-8") as output_handle:
        output_handle.write("##fileformat=VCFv4.3\n")
        output_handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        for record in sortable_records:
            _ensure_deletion_end(record)
            info_fields = record.get("INFO", {})
            serialized_info = _format_info(info_fields) if isinstance(info_fields, dict) else "."
            record_id = record.get("ID", ".")
            if isinstance(record_id, str):
                record_id = record_id.replace("%", "%25")
            output_handle.write(
                f"{record['CHROM']}\t{record['POS']}\t{record_id}\t{record['REF']}\t"
                f"{record['ALT']}\t.\tPASS\t{serialized_info}\n"
            )


def _resolve_vcf_midsv_tags(
    consensus_key: str,
    consensus_midsv_tags: list[str],
    control_sequence: str | None,
) -> list[str]:
    """Choose the MIDSV tokens that should be used for VCF export."""
    allele_name = _extract_allele_name(consensus_key)
    if not _should_align_to_control(allele_name) or not control_sequence:
        return consensus_midsv_tags

    consensus_sequence = convert_midsvs_to_sequence(consensus_midsv_tags).upper()
    aligned_midsv_tags = _align_sequence_to_control(control_sequence, consensus_sequence)
    return aligned_midsv_tags or consensus_midsv_tags


def export_to_vcf(
    tempdir: Path,
    sample_name: str,
    genome_coordinates: dict | None,
    cons_midsv_tags: dict[str, list[str]],
) -> None:
    """Export consensus MIDSV calls as per-allele VCF files."""
    normalized_genome_coordinates = genome_coordinates or {}
    chromosome_name = str(normalized_genome_coordinates.get("chrom") or "control")
    control_sequence = _load_control_sequence(tempdir, sample_name)

    for consensus_key, consensus_midsv_tags in cons_midsv_tags.items():
        record_id = consensus_key.replace("|", "_")
        midsv_tags_for_vcf = _resolve_vcf_midsv_tags(consensus_key, consensus_midsv_tags, control_sequence)
        local_records = _midsv_to_vcf_records(
            midsv_tags_for_vcf,
            chromosome_name,
            0,
            record_id,
            LARGE_SV_THRESHOLD,
        )
        genomic_records = _convert_records_to_genomic_coordinates(local_records, normalized_genome_coordinates)
        output_filename = build_report_filename(consensus_key, ".vcf", sample_name=sample_name)
        path_output = Path(tempdir, "report", "VCF", sample_name, output_filename)
        _write_vcf(path_output, genomic_records)
