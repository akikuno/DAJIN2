from __future__ import annotations

from typing import TypeAlias

FastqRecord: TypeAlias = dict[str, str]


def extract_qname_from_fastq_identifier(identifier: str) -> str:
    qname = identifier.split()[0]
    return qname[1:] if qname.startswith("@") else qname


def split_fastq_records_by_qname(
    fastq_records: list[FastqRecord], qnames_without_error: set[str]
) -> tuple[list[FastqRecord], list[FastqRecord]]:
    fastq_passed = []
    fastq_error = []
    for fastq_record in fastq_records:
        qname = extract_qname_from_fastq_identifier(fastq_record["identifier"])
        if qname in qnames_without_error:
            fastq_passed.append(fastq_record)
        else:
            fastq_error.append(fastq_record)

    return fastq_passed, fastq_error


def format_fastq_record(read: FastqRecord) -> str:
    return f"{read['identifier']}\n{read['sequence']}\n{read['separator']}\n{read['quality']}\n"
