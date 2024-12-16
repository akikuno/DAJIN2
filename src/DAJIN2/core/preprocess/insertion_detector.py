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


def convert_sv_tag_insertion(cssplits: list[str]) -> str:
    """Create a tag with I for bases that reads with insertion and M for bases that are mapped"""
    im_tags = []
    for tag in cssplits.split(","):
        if tag.startswith("+") and tag.count("+") > 10:
            ins = "I" * tag.count("|")
            im_tags.append(ins + "M")
        else:
            im_tags.append("M")

    return "".join(im_tags)


# TODO: def convert_sv_tag_deletion(cssplits: list[str]) -> str:
# TODO: def convert_sv_tag_inversion(cssplits: list[str]) -> str:


def convert_sv_tag(cssplits: list[str], sv_type: str) -> str:
    if sv_type == "insertion":
        return convert_sv_tag_insertion(cssplits)
    # TODO


def extract_sv_features(sv_tags: list[str], sv_type: str) -> np.ndarray[np.float64]:
    if sv_type == "insertion":
        key = "I"
    elif sv_type == "deletion":
        key = "D"
    elif sv_type == "inversion":
        key = "V"

    pattern = re.compile(rf"{key}+")

    features = []
    for tag in sv_tags:
        d_total_count = tag.count(key)
        d_ratio = d_total_count / len(tag) if len(tag) > 0 else 0
        group_count = len(pattern.findall(tag))
        features.append([d_total_count, d_ratio, group_count])

    return np.array(features)


def detect_deletion_alleles(
    TEMPDIR: Path, SAMPLE_NAME: str, CONTROL_NAME: str, FASTA_ALLELES: dict, sv_type: str
) -> None:
    path_midsv_sample = Path(TEMPDIR, SAMPLE_NAME, "midsv", "control", f"{SAMPLE_NAME}.jsonl")
    midsv_sample = io.read_jsonl(path_midsv_sample)
    sv_tags_sample = [convert_sv_tag(m["CSSPLIT"], sv_type) for m in midsv_sample]

    path_midsv_control = Path(TEMPDIR, CONTROL_NAME, "midsv", "control", f"{CONTROL_NAME}.jsonl")
    midsv_control = io.read_jsonl(path_midsv_control)
    sv_tags_control = [convert_sv_tag(m["CSSPLIT"], sv_type) for m in midsv_control][:1000]

    sv_features_sample = extract_sv_features(sv_tags_sample, sv_type)
    sv_features_control = extract_sv_features(sv_tags_control, sv_type)

    X = np.concatenate([sv_features_sample, sv_features_control])
    coverage_control = len(sv_features_control)
    coverage_sample = len(sv_features_sample)

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

    cstag_by_label = add_unique_allele_keys(cstag_by_label, FASTA_ALLELES, key=sv_type)
    fasta_by_label = add_unique_allele_keys(fasta_by_label, FASTA_ALLELES, key=sv_type)

    save_cstag(TEMPDIR, SAMPLE_NAME, cstag_by_label)
    save_fasta(TEMPDIR, SAMPLE_NAME, fasta_by_label)
