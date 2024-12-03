from __future__ import annotations

import gzip
import random
from pathlib import Path

import numpy as np
from rapidfuzz.distance import JaroWinkler
from rapidfuzz.process import cdist
from scipy.sparse import hstack
from sklearn.cluster import KMeans
from sklearn.feature_extraction.text import TfidfVectorizer

from DAJIN2.utils import io
from DAJIN2.utils.fastx_handler import read_fastq

###############################################################################
# Detect sequence errors
###############################################################################


def convert_nm_tag(csv_tags: list[list[str]]) -> str:
    """Create a tag with N for bases that reads are truncated and M for bases that are mapped"""
    return "".join(["N" if tag.startswith("N") or tag.startswith("n") else "M" for tag in csv_tags])


def detect_sequence_error_reads_in_control(ARGS) -> None:
    # Convert CSV strings to MIDSV tags
    midsv_control = io.read_jsonl(
        Path(ARGS.tempdir, ARGS.control_name, "midsv", "control", f"{ARGS.control_name}.jsonl")
    )
    nm_tags, qnames = zip(*[(convert_nm_tag(m["CSSPLIT"].split(",")), m["QNAME"]) for m in midsv_control])

    # Vectorize the MIDSV tags using TF-IDF with character-level N-grams
    vectorizer = TfidfVectorizer(analyzer="char", ngram_range=(100, 100))
    X = vectorizer.fit_transform(nm_tags)
    # Add a feature for the number of matches in the X
    match_counts = np.array([tag.count("M") for tag in nm_tags], dtype=int)
    X = hstack([X, match_counts.reshape(-1, 1)])

    # Apply KMeans clustering for binary classification based on the similarity of MIDSV.
    kmeans = KMeans(n_clusters=2, random_state=1)
    labels = kmeans.fit_predict(X)

    # Initialize lists for match counts of each label
    match_counts = {0: [], 1: []}

    # Count occurrences of "M" for each label
    for midsv_tag, label in zip(nm_tags, labels):
        match_counts[label].append(midsv_tag.count("M"))

    # Determine which label corresponds to the sequences with fewer matches and treat them as sequence errors.
    error_label = 0 if np.median(match_counts[0]) < np.median(match_counts[1]) else 1

    # Collect QNAMEs of sequences with or without sequence errors
    qnames_with_sequence_error = [qname for qname, label in zip(qnames, labels) if label == error_label]
    qnames_without_sequence_error = [qname for qname, label in zip(qnames, labels) if label != error_label]

    error_fraction = len(qnames_with_sequence_error) / len(qnames)

    # Output QNAMEs of sequences with sequence errors and the fraction of sequence errors
    Path(ARGS.tempdir, ARGS.control_name, "sequence_error").mkdir(parents=True, exist_ok=True)
    path_qnames_with_sequence_error = Path(
        ARGS.tempdir, ARGS.control_name, "sequence_error", "qnames_with_sequence_error.txt"
    )
    path_qnames_with_sequence_error.write_text("\n".join(qnames_with_sequence_error) + "\n")
    path_qnames_without_sequence_error = Path(
        ARGS.tempdir, ARGS.control_name, "sequence_error", "qnames_without_sequence_error.txt"
    )
    path_qnames_without_sequence_error.write_text("\n".join(qnames_without_sequence_error) + "\n")
    path_error_fraction = Path(ARGS.tempdir, ARGS.control_name, "sequence_error", "error_fraction.txt")
    path_error_fraction.write_text(str(error_fraction) + "\n")


def detect_sequence_error_reads_in_sample(ARGS) -> None:
    path_qnames_without_sequence_error = Path(
        ARGS.tempdir, ARGS.control_name, "sequence_error", "qnames_with_sequence_error.txt"
    )
    qnames_with_sequence_error_control = set(path_qnames_without_sequence_error.read_text().splitlines())

    midsv_control = io.read_jsonl(
        Path(ARGS.tempdir, ARGS.control_name, "midsv", "control", f"{ARGS.control_name}.jsonl")
    )
    midsv_errors = (m for m in midsv_control if m["QNAME"] in qnames_with_sequence_error_control)
    nm_tags_error = [convert_nm_tag(m["CSSPLIT"].split(",")) for m in midsv_errors]

    # ランダムに100本のエラー配列を取得
    random.seed(1)
    nm_tags_error = random.sample(nm_tags_error, min(len(nm_tags_error), 100))

    path_midsv_sample = Path(ARGS.tempdir, ARGS.sample_name, "midsv", "control", f"{ARGS.sample_name}.jsonl")
    midsv_sample = io.read_jsonl(path_midsv_sample)
    nm_tags_sample = [convert_nm_tag(m["CSSPLIT"].split(",")) for m in midsv_sample]

    similarity_scores = cdist(nm_tags_sample, nm_tags_error, scorer=JaroWinkler.normalized_similarity)
    most_similar_scores = np.max(similarity_scores, axis=1)

    midsv_sample = io.read_jsonl(path_midsv_sample)
    qnames_without_sequence_error = [m["QNAME"] for m, score in zip(midsv_sample, most_similar_scores) if score < 0.99]
    # Output QNAMEs of sequences with sequence errors and the fraction of sequence errors
    Path(ARGS.tempdir, ARGS.sample_name, "sequence_error").mkdir(parents=True, exist_ok=True)
    path_qnames_without_sequence_error = Path(
        ARGS.tempdir, ARGS.sample_name, "sequence_error", "qnames_without_sequence_error.txt"
    )
    path_qnames_without_sequence_error.write_text("\n".join(qnames_without_sequence_error) + "\n")


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

    path_qnames_without_sequence_error = Path(
        ARGS.tempdir, NAME, "sequence_error", "qnames_without_sequence_error.txt"
    )
    qnames_without_error = set(path_qnames_without_sequence_error.read_text().splitlines())

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
    path_qnames_without_sequence_error = Path(
        ARGS.tempdir, ARGS.sample_name, "sequence_error", "qnames_without_sequence_error.txt"
    )
    qnames_without_sequence_error = set(path_qnames_without_sequence_error.read_text().splitlines())

    for path_midsv in Path(ARGS.tempdir, ARGS.sample_name, "midsv").glob(f"*/{ARGS.sample_name}.jsonl"):
        midsv_sample = io.read_jsonl(path_midsv)
        midsv_sample_filtered = [m for m in midsv_sample if m["QNAME"] in qnames_without_sequence_error]
        io.write_jsonl(midsv_sample_filtered, path_midsv)
