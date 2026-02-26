from __future__ import annotations

import hashlib

MAX_FILENAME_BYTES = 255


def _truncate_utf8(value: str, max_bytes: int) -> str:
    """Trim a string to a UTF-8 byte limit without breaking multibyte characters."""
    encoded = value.encode("utf-8")
    if len(encoded) <= max_bytes:
        return value
    trimmed = encoded[:max_bytes]
    while trimmed and (trimmed[-1] & 0xC0) == 0x80:
        trimmed = trimmed[:-1]
    return trimmed.decode("utf-8", errors="ignore")


def _compact_consensus_label(label: str) -> str:
    """Shorten consensus labels by keeping allele id and percentage when possible."""
    normalized = label.replace("|", "_")
    parts = normalized.split("_")
    if len(parts) >= 2 and parts[0].startswith("allele") and parts[-1].endswith("%"):
        return f"{parts[0]}_{parts[-1]}"
    return normalized


def build_report_filename(
    label: str,
    extension: str,
    sample_name: str | None = None,
    max_filename_bytes: int = MAX_FILENAME_BYTES,
) -> str:
    """
    Build a safe report filename with compaction and hash fallback for long names to
    avoid 'OSError: [Errno 36] File name too long'.
    Example: "sample1_allele01_{long_allele_name}_12.3%.fasta" -> "sample1_allele01_12.3%.fasta" (if under limit).
    """
    normalized = label.replace("|", "_")
    prefix = f"{sample_name}_" if sample_name else ""
    filename = f"{prefix}{normalized}{extension}"
    if len(filename.encode("utf-8")) <= max_filename_bytes:
        return filename

    compact = _compact_consensus_label(label)
    filename = f"{prefix}{compact}{extension}"
    if len(filename.encode("utf-8")) <= max_filename_bytes:
        return filename

    digest = hashlib.md5(normalized.encode("utf-8")).hexdigest()[:10]
    suffix = f"_{digest}{extension}"
    max_stem_bytes = max_filename_bytes - len(suffix.encode("utf-8"))
    stem = _truncate_utf8(f"{prefix}{compact}", max_stem_bytes)
    return f"{stem}{suffix}"
