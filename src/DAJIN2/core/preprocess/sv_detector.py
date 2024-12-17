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


###############################################################################
# extract_sv_features
###############################################################################


def _find_sv_start_index(sv_tag: str, sv_type: str) -> list[int]:
    if sv_type == "insertion":
        key = "I"
    elif sv_type == "deletion":
        key = "D"
    elif sv_type == "inversion":
        key = "V"

    pattern = re.compile(rf"{key}+")
    return [match.start() for match in pattern.finditer(sv_tag)]


def _extract_sv_start_index(sv_tags: list[str], sv_type: str) -> list[int]:
    sv_start_index = []
    for tag in sv_tags:
        for i in _find_sv_start_index(tag, sv_type):
            sv_start_index.append(i)
    return sorted(sv_start_index)


def _group_sv_start_index(sv_tags: list[list[str]], sv_type: str) -> list[list[int]]:
    sv_start_index = _extract_sv_start_index(sv_tags, sv_type)

    groups = []
    current_group = []
    for key in sv_start_index:
        if not current_group or key - current_group[-1] <= 5:
            current_group.append(key)
        else:
            groups.append(current_group)
            current_group = [key]

    if current_group:
        groups.append(current_group)

    return [group for group in groups if len(group) > 10]


def define_index_converter(sv_tags: list[list[str]], sv_type: str) -> dict[int, int]:
    index_group = _group_sv_start_index(sv_tags, sv_type)
    index_converter = {}
    for group in index_group:
        value = group[0]
        index_converter |= {g: value for g in group}
    return index_converter


def annotate_merged_sv_start_index(
    sv_tags: list[str], sv_type: str, index_converter: dict[int, int]
) -> list[list[int]]:
    sv_start_index = []
    for tag in sv_tags:
        start_index = []
        for i in _find_sv_start_index(tag, sv_type):
            start_index.append(index_converter.get(i, None))
        sv_start_index.append(start_index)
    return sv_start_index


def extract_sv_start_index(sv_tags_index: list[list[int]], sv_index: list[int]) -> list[list[int]]:
    X_index = []
    for index in sv_tags_index:
        x_index = [0] * len(sv_index)
        for i in index:
            if i:
                x_index[sv_index.index(i)] = 1
        X_index.append(x_index)
    return X_index


def _find_sv_start_index_and_size(sv_tag: str, sv_type: str) -> list[int]:
    if sv_type == "insertion":
        key = "I"
    elif sv_type == "deletion":
        key = "D"
    elif sv_type == "inversion":
        key = "V"

    pattern = re.compile(rf"{key}+")
    return [[match.start(), match.end() - match.start()] for match in pattern.finditer(sv_tag)]


def annotate_size_by_merged_sv_start_index(
    sv_tags: list[str], sv_type: str, index_converter: dict[int, int]
) -> list[list[int]]:
    sv_start_index_size = []
    for tag in sv_tags:
        start_index_size = []
        for i, size in _find_sv_start_index_and_size(tag, sv_type):
            start_index_size.append([index_converter.get(i, None), size])
        sv_start_index_size.append(start_index_size)
    return sv_start_index_size


def extract_sv_size(sv_tags_size: list[list[list[int]]], sv_index: list[int]) -> list[list[int]]:
    X_size = []
    for index_size in sv_tags_size:
        x_size = [0] * len(sv_index)
        for i, size in index_size:
            if i:
                x_size[sv_index.index(i)] = size
        X_size.append(x_size)
    return X_size


def extract_sv_features(sv_tags: list[str], sv_type: str, sv_index: list[int]) -> np.ndarray[np.float64]:
    index_converter = define_index_converter(sv_tags, sv_type)
    sv_tags_index = annotate_merged_sv_start_index(sv_tags, sv_type, index_converter)
    sv_tags_size = annotate_size_by_merged_sv_start_index(sv_tags, sv_type, index_converter)

    X_index = extract_sv_start_index(sv_tags_index, sv_index)
    X_size = extract_sv_size(sv_tags_size, sv_index)

    return np.concatenate([X_index, X_size], axis=1)

###############################################################################
# 変異部についてコンセンサス配列を入手する
###############################################################################


###############################################################################
# main
###############################################################################

def detect_sv_alleles(TEMPDIR: Path, SAMPLE_NAME: str, CONTROL_NAME: str, FASTA_ALLELES: dict, sv_type: str) -> None:
    path_midsv_sample = Path(TEMPDIR, SAMPLE_NAME, "midsv", "control", f"{SAMPLE_NAME}.jsonl")
    midsv_sample = io.read_jsonl(path_midsv_sample)
    sv_tags_sample = [convert_sv_tag(m["CSSPLIT"], sv_type) for m in midsv_sample]

    path_midsv_control = Path(TEMPDIR, CONTROL_NAME, "midsv", "control", f"{CONTROL_NAME}.jsonl")
    midsv_control = io.read_jsonl(path_midsv_control)
    sv_tags_control = [convert_sv_tag(m["CSSPLIT"], sv_type) for m in midsv_control][:1000]

    sv_index = sorted(set(define_index_converter(sv_tags_sample, sv_type).values()))

    sv_features_sample = extract_sv_features(sv_tags_sample, sv_type, sv_index)
    sv_features_control = extract_sv_features(sv_tags_control, sv_type, sv_index)

    X = np.concatenate([sv_features_sample, sv_features_control])
    coverage_control = len(sv_features_control)
    coverage_sample = len(sv_features_sample)

    labels = optimize_labels(X, coverage_sample, coverage_control)

    midsv_sample = io.read_jsonl(path_midsv_sample)
    cssplits_iter: Iterator[list[list[str]]] = (m["CSSPLIT"].split(",") for m in midsv_sample)

    cssplits_by_label = defaultdict(list)
    for label, cssplits in zip(labels, cssplits_iter):
        cssplits_by_label[label].append(cssplits)

    #######################################################
    # Generate cstag consensus
    #######################################################
    cstag_by_label = {}
    fasta_by_label = {}
    for label, cssplits in cssplits_by_label.items():
        cssplits_subset = cssplits[:100]
        cstag_subset = [convert_cssplits_to_cstag(tag) for tag in cssplits_subset]
        positions = [1] * len(cstag_subset)
        cstag_consensus = cstag.consensus(cstag_subset, positions)
        cstag_by_label[label] = cstag_consensus
        fasta_by_label[label] = cstag.to_sequence(cstag_consensus)

    #######################################################
    # Remove similar alleles to user's alleles, or clustered alleles
    #######################################################
    fasta_by_label = extract_unique_sv(fasta_by_label, FASTA_ALLELES)
    cstag_by_label = {label: cs_tag for label, cs_tag in cstag_by_label.items() if label in fasta_by_label}

    #######################################################
    # Output cstag and fasta
    #######################################################
    cstag_by_label = add_unique_allele_keys(cstag_by_label, FASTA_ALLELES, key=sv_type)
    fasta_by_label = add_unique_allele_keys(fasta_by_label, FASTA_ALLELES, key=sv_type)

    save_cstag(TEMPDIR, SAMPLE_NAME, cstag_by_label)
    save_fasta(TEMPDIR, SAMPLE_NAME, fasta_by_label)
