from __future__ import annotations

from pathlib import Path


def create_temporal(TEMPDIR: Path, NAME: str, is_control=False) -> None:
    Path(TEMPDIR, "result").mkdir(parents=True, exist_ok=True)
    if is_control:
        SUBDIRS = ["fasta", "sam", "midsv", "mutation_loci", "clustering"]
    else:
        SUBDIRS = [
            "fasta",
            "sam",
            "midsv",
            "mutation_loci",
            "knockin_loci",
            "classification",
            "clustering",
            "consensus",
        ]
    for subdir in SUBDIRS:
        Path(TEMPDIR, NAME, subdir).mkdir(parents=True, exist_ok=True)


def create_report(TEMPDIR: Path, NAME: str, is_control=False) -> None:
    if is_control:
        Path(TEMPDIR, "report", "BAM", NAME).mkdir(parents=True, exist_ok=True)
        return
    SUBDIRS_REPORT = ["HTML", "FASTA", "BAM", "MUTATION_INFO", ".igvjs"]
    for reportdir in SUBDIRS_REPORT:
        if reportdir == "MUTATION_INFO":
            Path(TEMPDIR, "report", reportdir).mkdir(parents=True, exist_ok=True)
        else:
            Path(TEMPDIR, "report", reportdir, NAME).mkdir(parents=True, exist_ok=True)
