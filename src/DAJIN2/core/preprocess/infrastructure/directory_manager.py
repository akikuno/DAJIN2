from __future__ import annotations

from pathlib import Path


def create_temporal_directories(TEMPDIR: Path, NAME: str, is_control=False) -> None:
    Path(TEMPDIR, "result").mkdir(parents=True, exist_ok=True)
    SUBDIRS = ["fasta", "fastq", "sam", "midsv", "mutation_loci", "clustering", "consensus"]
    if is_control is False:
        SUBDIRS.extend(["cstag", "knockin_loci", "classification"])
    for subdir in SUBDIRS:
        Path(TEMPDIR, NAME, subdir).mkdir(parents=True, exist_ok=True)


def create_report_directories(TEMPDIR: Path, NAME: str, is_control=False) -> None:
    if is_control:
        Path(TEMPDIR, "report", "BAM", NAME).mkdir(parents=True, exist_ok=True)
        return
    SUBDIRS_REPORT = ["HTML", "FASTA", "BAM", "MUTATION_INFO", ".igvjs"]
    for reportdir in SUBDIRS_REPORT:
        if reportdir == "MUTATION_INFO":
            Path(TEMPDIR, "report", reportdir).mkdir(parents=True, exist_ok=True)
        else:
            Path(TEMPDIR, "report", reportdir, NAME).mkdir(parents=True, exist_ok=True)
