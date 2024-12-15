from __future__ import annotations

import re
from collections import defaultdict
from collections.abc import Iterator
from pathlib import Path

import cstag
import numpy as np

from DAJIN2.core.clustering.clustering import optimize_labels
from DAJIN2.core.preprocess.sv_handler import add_unique_allele_keys, extract_unique_sv, save_cstag, save_fasta
from DAJIN2.utils import io
from DAJIN2.utils.cssplits_handler import convert_cssplits_to_cstag

###############################################################################
# Detect sequence errors
###############################################################################


def convert_dm_tag(csv_tags: list[list[str]]) -> str:
    """Create a tag with D for bases that reads are truncated and M for bases that are mapped"""
    return "".join(["D" if tag.startswith("-") else "M" for tag in csv_tags])


d_pattern = re.compile(r"D+")


def extract_d_features(dm_tags: list[str]) -> np.ndarray[np.float64]:
    features = []
    for tag in dm_tags:
        d_total_count = tag.count("D")
        d_ratio = d_total_count / len(tag) if len(tag) > 0 else 0
        group_count = len(d_pattern.findall(tag))
        features.append([d_ratio, group_count])
    return np.array(features)


def detect_deletion_alleles(TEMPDIR: Path, SAMPLE_NAME: str, CONTROL_NAME: str, FASTA_ALLELES: dict) -> None:
    path_midsv_sample = Path(TEMPDIR, SAMPLE_NAME, "midsv", "control", f"{SAMPLE_NAME}.jsonl")
    midsv_sample = io.read_jsonl(path_midsv_sample)
    dm_tags_sample = [convert_dm_tag(m["CSSPLIT"].split(",")) for m in midsv_sample]

    path_midsv_control = Path(TEMPDIR, CONTROL_NAME, "midsv", "control", f"{CONTROL_NAME}.jsonl")
    midsv_control = io.read_jsonl(path_midsv_control)
    dm_tags_control = [convert_dm_tag(m["CSSPLIT"].split(",")) for m in midsv_control][:1000]

    d_features_sample = extract_d_features(dm_tags_sample)
    d_features_control = extract_d_features(dm_tags_control)

    X = np.concatenate([d_features_sample, d_features_control])
    coverage_control = len(d_features_control)
    coverage_sample = len(d_features_sample)

    labels = optimize_labels(X, coverage_sample, coverage_control)

    midsv_sample = io.read_jsonl(path_midsv_sample)
    cssplits_iter: Iterator[list[list[str]]] = (m["CSSPLIT"].split(",") for m in midsv_sample)

    cssplits_by_label = defaultdict(list)
    for label, cssplits in zip(labels, cssplits_iter):
        cssplits_by_label[label].append(cssplits)

    cstag_by_label = {}
    fasta_by_label = {}
    for label, cssplits in cssplits_by_label.items():
        cssplits_subset = cssplits[:100]
        cstag_subset = [convert_cssplits_to_cstag(tag) for tag in cssplits_subset]
        positions = [1] * len(cstag_subset)
        cstag_consensus = cstag.consensus(cstag_subset, positions)
        cstag_by_label[label] = cstag_consensus
        fasta_by_label[label] = cstag.to_sequence(cstag_consensus)

    # Remove similar alleles to user's alleles, or clustered alleles
    fasta_by_label = extract_unique_sv(fasta_by_label, FASTA_ALLELES)
    cstag_by_label = {label: cs_tag for label, cs_tag in cstag_by_label.items() if label in fasta_by_label}

    cstag_by_label = add_unique_allele_keys(cstag_by_label, FASTA_ALLELES, key="deletion")
    fasta_by_label = add_unique_allele_keys(fasta_by_label, FASTA_ALLELES, key="deletion")

    save_cstag(TEMPDIR, SAMPLE_NAME, cstag_by_label)
    save_fasta(TEMPDIR, SAMPLE_NAME, fasta_by_label)
