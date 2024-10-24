from __future__ import annotations

import gzip
from collections import Counter
from pathlib import Path

import numpy as np
from sklearn.cluster import HDBSCAN, KMeans
from sklearn.feature_extraction.text import TfidfVectorizer

from DAJIN2.utils import io
from DAJIN2.utils.fastx_handler import parse_fastq, read_lines

###############################################################################
# Detect sequence errors
###############################################################################


def parse_midsv_from_csv(csv_tags: list[list[str]]) -> str:
    midsv_seq = []
    for tag in csv_tags:
        if tag.startswith("="):
            midsv_seq.append("M")
        elif tag.startswith("+"):
            midsv_seq.append("I")
        elif tag.startswith("-"):
            midsv_seq.append("D")
        elif tag.startswith("*"):
            midsv_seq.append("S")
        elif tag.startswith("N"):
            midsv_seq.append("N")
    return "".join(midsv_seq)


def detect_sequence_error_reads_in_control(ARGS) -> None:
    # Convert CSV strings to MIDSV tags
    midsv_control = io.read_jsonl(Path(ARGS.tempdir, ARGS.control_name, "midsv", "control", "control.jsonl"))
    midsv_tags, qnames = zip(*[(parse_midsv_from_csv(m["CSSPLIT"].split(",")), m["QNAME"]) for m in midsv_control])

    # Vectorize the MIDSV tags using TF-IDF with character-level 3-grams
    vectorizer = TfidfVectorizer(analyzer="char", ngram_range=(3, 3))
    X = vectorizer.fit_transform(midsv_tags)

    # Apply KMeans clustering for binary classification based on the similarity of MIDSV.
    kmeans = KMeans(n_clusters=2, random_state=1)
    labels = kmeans.fit_predict(X)

    # Initialize lists for match counts of each label
    match_counts = {0: [], 1: []}

    # Count occurrences of "M" for each label
    for midsv_tag, label in zip(midsv_tags, labels):
        match_counts[label].append(midsv_tag.count("M"))

    # Determine which label corresponds to the sequences with fewer matches and treat them as sequence errors.
    error_label = 0 if np.median(match_counts[0]) < np.median(match_counts[1]) else 1

    # Collect QNAMEs of sequences with sequence errors
    qnames_sequence_error = [qname for qname, label in zip(qnames, labels) if label == error_label]

    error_fraction = len(qnames_sequence_error) / len(qnames)

    # Output QNAMEs of sequences with sequence errors and the fraction of sequence errors
    Path(ARGS.tempdir, ARGS.control_name, "sequence_error").mkdir(parents=True, exist_ok=True)
    path_qnames_error_control = Path(ARGS.tempdir, ARGS.control_name, "sequence_error", "qnames_sequence_error.txt")
    path_qnames_error_control.write_text("\n".join(qnames_sequence_error) + "\n")
    path_error_fraction = Path(ARGS.tempdir, ARGS.control_name, "sequence_error", "error_fraction.txt")
    path_error_fraction.write_text(str(error_fraction) + "\n")


def detect_sequence_error_reads_in_sample(ARGS) -> None:
    path_qnames_error_control = Path(ARGS.tempdir, ARGS.control_name, "sequence_error", "qnames_sequence_error.txt")
    qnames_error_control = set(path_qnames_error_control.read_text().splitlines())
    path_error_fraction_control = Path(ARGS.tempdir, ARGS.control_name, "sequence_error", "error_fraction.txt")
    error_fraction_control = float(path_error_fraction_control.read_text())

    midsv_control = io.read_jsonl(Path(ARGS.tempdir, ARGS.control_name, "midsv", "control", "control.jsonl"))
    midsv_errors = (m for m in midsv_control if m["QNAME"] in qnames_error_control)
    midsv_tags_error = [parse_midsv_from_csv(m["CSSPLIT"].split(",")) for m in midsv_errors]

    path_midsv_sample = Path(ARGS.tempdir, ARGS.sample_name, "midsv", "control", f"{ARGS.sample_name}.jsonl")
    midsv_sample = io.read_jsonl(path_midsv_sample)
    midsv_tags_sample = [parse_midsv_from_csv(m["CSSPLIT"].split(",")) for m in midsv_sample]

    midsv_tags_all = midsv_tags_error + midsv_tags_sample

    sample_number = io.count_newlines(path_midsv_sample)
    min_cluster_size = int(sample_number * error_fraction_control)

    vectorizer = TfidfVectorizer(analyzer="char", ngram_range=(3, 3)).fit(midsv_tags_all)
    X_all = vectorizer.transform(midsv_tags_all)
    labels_hdbscan = HDBSCAN(min_cluster_size=min_cluster_size).fit_predict(X_all).tolist()
    len_seq_error = len(midsv_tags_error)
    labels_error = labels_hdbscan[:len_seq_error]
    labels_sample = labels_hdbscan[len_seq_error:]
    label_error = Counter(labels_error).most_common()[0][0]

    midsv_sample = io.read_jsonl(path_midsv_sample)
    error_qnames_sample = [m["QNAME"] for m, label in zip(midsv_sample, labels_sample) if label == label_error]

    # Output QNAMEs of sequences with sequence errors and the fraction of sequence errors
    Path(ARGS.tempdir, ARGS.sample_name, "sequence_error").mkdir(parents=True, exist_ok=True)
    path_qnames_error_control = Path(ARGS.tempdir, ARGS.sample_name, "sequence_error", "qnames_sequence_error.txt")
    path_qnames_error_control.write_text("\n".join(error_qnames_sample) + "\n")


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

    path_qnames_error_sample = Path(ARGS.tempdir, NAME, "sequence_error", "qnames_sequence_error.txt")
    path_fastq = Path(ARGS.tempdir, NAME, "fastq", f"{NAME}.fastq.gz")
    path_fastq_error = Path(ARGS.tempdir, NAME, "fastq", f"{NAME}_sequence_error.fastq.gz")
    qnames_error_sample = set(path_qnames_error_sample.read_text().splitlines())

    fastq: list[dict] = parse_fastq(read_lines(path_fastq))

    # -----------------------------------------------------
    # Split FASTQ by sequence error
    # -----------------------------------------------------
    fastq_filtered = []
    fastq_error = []
    for fastq_record in fastq:
        qname = fastq_record["header"].split()[0][1:]
        if qname not in qnames_error_sample:
            fastq_filtered.append(fastq_record)
        else:
            fastq_error.append(fastq_record)

    # -----------------------------------------------------
    # Output FASTQ files
    # -----------------------------------------------------
    with gzip.open(path_fastq, "wt") as f:
        for read in fastq_filtered:
            f.write(f"{read['header']}\n{read['sequence']}\n{read['annotate']}\n{read['quality']}\n")

    with gzip.open(path_fastq_error, "wt") as f:
        for read in fastq_error:
            f.write(f"{read['header']}\n{read['sequence']}\n{read['annotate']}\n{read['quality']}\n")


###############################################################################
# Replace MIDSV files without sequence errors
###############################################################################


def replace_midsv_without_sequence_errors(ARGS) -> None:
    path_qnames_error_sample = Path(ARGS.tempdir, ARGS.sample_name, "sequence_error", "qnames_sequence_error.txt")
    qnames_error_sample = set(path_qnames_error_sample.read_text().splitlines())

    for path_midsv in Path(ARGS.tempdir, ARGS.sample_name, "midsv").glob(f"*/{ARGS.sample_name}.jsonl"):
        midsv_sample = io.read_jsonl(path_midsv)
        midsv_sample_filtered = [m for m in midsv_sample if m["QNAME"] not in qnames_error_sample]
        io.write_jsonl(midsv_sample_filtered, path_midsv)
