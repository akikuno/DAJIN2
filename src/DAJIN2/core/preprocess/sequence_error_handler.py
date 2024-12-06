from __future__ import annotations

import gzip
from pathlib import Path
from typing import Callable

import numpy as np
from sklearn.cluster import KMeans
from sklearn.linear_model import LogisticRegression

from DAJIN2.utils import io
from DAJIN2.utils.fastx_handler import read_fastq

###############################################################################
# Detect sequence errors
###############################################################################


def convert_nm_tag(csv_tags: list[list[str]]) -> str:
    """Create a tag with N for bases that reads are truncated and M for bases that are mapped"""
    return "".join(["N" if tag.startswith("N") or tag.startswith("n") else "M" for tag in csv_tags])


def extract_n_features(nm_tags: list[str]) -> np.ndarray[np.float64]:
    features = []
    for s in nm_tags:
        total_length = len(s)
        n_count = s.count("N")
        n_ratio = n_count / total_length if total_length > 0 else 0
        max_n_run = max([len(run) for run in s.split("M")] + [0])
        features.append([n_ratio, max_n_run])
    return np.array(features)


def detect_sequence_error_reads_in_control(ARGS) -> None:
    # Convert CSV strings to MIDSV tags
    midsv_control = io.read_jsonl(
        Path(ARGS.tempdir, ARGS.control_name, "midsv", "control", f"{ARGS.control_name}.jsonl")
    )
    nm_tags, qnames = zip(*[(convert_nm_tag(m["CSSPLIT"].split(",")), m["QNAME"]) for m in midsv_control])

    X = extract_n_features(nm_tags)
    kmeans = KMeans(n_clusters=2, random_state=1)
    labels = kmeans.fit_predict(X)

    # Determine which label corresponds to the sequences with fewer matches and treat them as sequence errors.
    match_counts = {0: [], 1: []}
    for nm_tag, label in zip(nm_tags, labels):
        match_counts[label].append(nm_tag.count("M"))

    error_label = 0 if np.median(match_counts[0]) < np.median(match_counts[1]) else 1

    # Collect QNAMEs of sequences with or without sequence errors
    qnames_with_error = [qname for qname, label in zip(qnames, labels) if label == error_label]
    qnames_without_error = [qname for qname, label in zip(qnames, labels) if label != error_label]

    error_fraction = len(qnames_with_error) / len(qnames)

    # Output QNAMEs of sequences with sequence errors and the fraction of sequence errors
    Path(ARGS.tempdir, ARGS.control_name, "sequence_error").mkdir(parents=True, exist_ok=True)

    # Output QNAMEs of sequences with sequence errors and the fraction of sequence errors
    path_qnames_without_error = Path(ARGS.tempdir, ARGS.control_name, "sequence_error", "qnames_without_error.txt")
    path_qnames_without_error.write_text("\n".join(qnames_without_error) + "\n")
    path_error_fraction = Path(ARGS.tempdir, ARGS.control_name, "sequence_error", "error_fraction.txt")
    path_error_fraction.write_text(str(error_fraction) + "\n")


def load_midsv_to_nm_tags(file_path: Path, filter_condition: Callable) -> list[str]:
    """Load MIDSV file and convert to NM tags"""
    midsv = io.read_jsonl(file_path)
    filtered_midsv = (m for m in midsv if filter_condition(m))
    return [convert_nm_tag(m["CSSPLIT"].split(",")) for m in filtered_midsv]


def detect_sequence_error_reads_in_sample(ARGS) -> None:
    path_qnames_without_error = Path(ARGS.tempdir, ARGS.control_name, "sequence_error", "qnames_without_error.txt")
    qnames_without_error = set(path_qnames_without_error.read_text().splitlines())

    path_midsv_control = Path(ARGS.tempdir, ARGS.control_name, "midsv", "control", f"{ARGS.control_name}.jsonl")
    nm_tags_without_error = load_midsv_to_nm_tags(path_midsv_control, lambda m: m["QNAME"] in qnames_without_error)
    nm_tags_with_error: list[str] = load_midsv_to_nm_tags(
        path_midsv_control,
        lambda m: m["QNAME"] not in qnames_without_error,
    )

    n_features_without_error = extract_n_features(nm_tags_without_error)
    n_features_with_error = extract_n_features(nm_tags_with_error)

    X = np.vstack([n_features_without_error, n_features_with_error])
    y = np.array([0] * len(n_features_without_error) + [1] * len(n_features_with_error))
    clf = LogisticRegression(random_state=0).fit(X, y)

    path_midsv_sample = Path(ARGS.tempdir, ARGS.sample_name, "midsv", "control", f"{ARGS.sample_name}.jsonl")
    midsv_sample = io.read_jsonl(path_midsv_sample)
    nm_tags_sample = [convert_nm_tag(m["CSSPLIT"].split(",")) for m in midsv_sample]

    n_features_sample = extract_n_features(nm_tags_sample)
    labels = clf.predict(n_features_sample)

    midsv_sample = io.read_jsonl(path_midsv_sample)
    qnames_without_error = [m["QNAME"] for m, label in zip(midsv_sample, labels) if label == 0]

    # Output QNAMEs of sequences with sequence errors and the fraction of sequence errors
    Path(ARGS.tempdir, ARGS.sample_name, "sequence_error").mkdir(parents=True, exist_ok=True)
    path_qnames_without_error = Path(ARGS.tempdir, ARGS.sample_name, "sequence_error", "qnames_without_error.txt")
    path_qnames_without_error.write_text("\n".join(qnames_without_error) + "\n")


def detect_sequence_error_reads(ARGS, is_control: bool = False) -> None:
    if is_control:
        detect_sequence_error_reads_in_control(ARGS)
    else:
        detect_sequence_error_reads_in_sample(ARGS)


###############################################################################
# Split FASTQ by sequence error
###############################################################################


def split_fastq_by_sequence_error(ARGS, is_control: bool = False) -> None:
    if is_control:
        NAME = ARGS.control_name
    else:
        NAME = ARGS.sample_name

    path_qnames_without_error = Path(ARGS.tempdir, NAME, "sequence_error", "qnames_without_error.txt")
    qnames_without_error = set(path_qnames_without_error.read_text().splitlines())

    path_fastq = Path(ARGS.tempdir, NAME, "fastq", f"{NAME}.fastq.gz")
    path_fastq_error = Path(ARGS.tempdir, NAME, "fastq", f"{NAME}_sequence_error.fastq.gz")

    fastq: list[dict] = read_fastq(path_fastq)

    # -----------------------------------------------------
    # Split FASTQ by sequence error
    # -----------------------------------------------------
    fastq_passed = []
    fastq_error = []
    for fastq_record in fastq:
        qname = fastq_record["identifier"].split()[0][1:]
        if qname in qnames_without_error:
            fastq_passed.append(fastq_record)
        else:
            fastq_error.append(fastq_record)

    # -----------------------------------------------------
    # Output FASTQ files
    # -----------------------------------------------------
    with gzip.open(path_fastq, "wt") as f:
        for read in fastq_passed:
            f.write(f"{read['identifier']}\n{read['sequence']}\n{read['separator']}\n{read['quality']}\n")

    with gzip.open(path_fastq_error, "wt") as f:
        for read in fastq_error:
            f.write(f"{read['identifier']}\n{read['sequence']}\n{read['separator']}\n{read['quality']}\n")


###############################################################################
# Replace MIDSV files without sequence errors
###############################################################################


def replace_midsv_without_sequence_errors(ARGS) -> None:
    path_qnames_without_error = Path(ARGS.tempdir, ARGS.sample_name, "sequence_error", "qnames_without_error.txt")
    qnames_without_error = set(path_qnames_without_error.read_text().splitlines())

    for path_midsv in Path(ARGS.tempdir, ARGS.sample_name, "midsv").glob(f"*/{ARGS.sample_name}.jsonl"):
        midsv_sample = io.read_jsonl(path_midsv)
        midsv_sample_filtered = [m for m in midsv_sample if m["QNAME"] in qnames_without_error]
        io.write_jsonl(midsv_sample_filtered, path_midsv)
