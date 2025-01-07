from __future__ import annotations

from collections import Counter, defaultdict
from collections.abc import Iterator
from itertools import chain
from pathlib import Path

import cstag
import numpy as np

from DAJIN2.core.clustering.clustering import optimize_labels
from DAJIN2.core.preprocess.sv_handler import add_unique_allele_keys, extract_unique_sv, save_cstag, save_fasta
from DAJIN2.utils import io
from DAJIN2.utils.cssplits_handler import convert_cssplits_to_cstag, revcomp_cssplits

###############################################################################
# Clusteringにつかう特徴量を抽出する
# 特徴量：SVの開始のIndexとSVの大きさ
###############################################################################


def get_index_of_sv(misdv_string: str, sv_type: str, mutation_loci: list[set[str]]) -> list[int]:
    """Return a list with the start index of SV"""

    items = misdv_string.split(",")
    index_of_sv = []
    previous_midsv = items[0]

    for i, item in enumerate(items):
        if sv_type in ("deletion", "inversion"):
            if sv_type == "deletion":
                is_sv_item = item.startswith("-") and "-" in mutation_loci[i]
                is_prev_sv_item = not previous_midsv.startswith("-")
            else:  # inversion
                is_sv_item = item.islower()
                is_prev_sv_item = not previous_midsv.islower()

            if is_sv_item and is_prev_sv_item:
                index_of_sv.append(i)

        elif sv_type == "insertion":
            if item.startswith("+") and "+" in mutation_loci[i]:
                index_of_sv.append(i)

        previous_midsv = item
    return index_of_sv


def group_similar_indices(index_of_sv: list[int], distance: int = 5, min_coverage: float = 5.0) -> list[list[int]]:
    """Return a list of groups of similar indices"""
    index_of_sv.sort()
    groups = []
    current_group = []
    for key in index_of_sv:
        if not current_group or key - current_group[-1] <= distance:
            current_group.append(key)
        else:
            groups.append(current_group)
            current_group = [key]

    if current_group:
        groups.append(current_group)

    return [group for group in groups if len(group) > min_coverage]


def define_index_converter(index_grouped: list[list[int]]) -> dict[int, int]:
    index_converter = {}
    for group in index_grouped:
        value = Counter(group).most_common()[0][0]
        index_converter |= {g: value for g in group}
    return index_converter


def get_index_and_sv_size(
    misdv_string: str, sv_type: str, mutation_loci: list[set[str]], index_converter: dict[int, int]
) -> dict[int, int]:
    """Return a dictionary with the start index and SV size"""

    items = misdv_string.split(",")
    result = {}
    sv_size = 0
    start_index = None
    previous_midsv = items[0]

    for i, item in enumerate(items):
        if sv_type in ("deletion", "inversion"):
            if sv_type == "deletion":
                is_sv_item = item.startswith("-") and "-" in mutation_loci[i]
                is_prev_sv_item = previous_midsv.startswith("-")
            else:  # inversion
                is_sv_item = item.islower()
                is_prev_sv_item = previous_midsv.islower()

            # Get the start index of SV
            if is_sv_item and not is_prev_sv_item and i in index_converter:
                start_index = index_converter[i]

            # Calculate the size of SV
            if is_sv_item:
                sv_size += 1
            elif start_index:
                result.update({start_index: sv_size})
                start_index = None
                sv_size = 0

        elif sv_type == "insertion":
            if item.startswith("+") and "+" in mutation_loci[i] and i in index_converter:
                result.update({index_converter[i]: item.count("|")})

        previous_midsv = item
    return result


def extract_features(index_and_sv_size: list[dict[int, int]], all_sv_index: set[str]) -> np.ndarray:
    """補正後のIndexにあるSVの大きさを特徴量として用いる"""
    sv_size_features = []

    for record in index_and_sv_size:
        sv_size = {i: 0 for i in all_sv_index}
        sv_size |= record
        sv_size_features.append(sv_size)

    index_of_sv = sorted({key for record in sv_size_features for key in record.keys()})
    return np.array([[record.get(i, 0) for i in index_of_sv] for record in sv_size_features])


###############################################################################
# main
###############################################################################


def detect_sv_alleles(TEMPDIR: Path, SAMPLE_NAME: str, CONTROL_NAME: str, FASTA_ALLELES: dict, sv_type: str) -> None:
    #######################################################
    # Load data
    #######################################################
    MUTATION_LOCI: list[set[str]] = io.load_pickle(
        Path(TEMPDIR, SAMPLE_NAME, "mutation_loci", "control", "mutation_loci.pickle")
    )

    path_midsv_sample = Path(TEMPDIR, SAMPLE_NAME, "midsv", "control", f"{SAMPLE_NAME}.jsonl")
    path_midsv_control = Path(TEMPDIR, CONTROL_NAME, "midsv", "control", f"{CONTROL_NAME}.jsonl")

    coverage_sample = io.count_newlines(path_midsv_sample)
    coverage_control = io.count_newlines(path_midsv_control)

    #######################################################
    # Extract features: SV start index and SV size
    #######################################################
    index_of_sv_sample = [
        get_index_of_sv(m["CSSPLIT"], sv_type, MUTATION_LOCI) for m in io.read_jsonl(path_midsv_sample)
    ]
    index_of_sv_sample = sorted(chain.from_iterable(index_of_sv_sample))

    index_of_sv_control = [
        get_index_of_sv(m["CSSPLIT"], sv_type, MUTATION_LOCI) for m in io.read_jsonl(path_midsv_control)
    ][:1000]
    index_of_sv_control = sorted(chain.from_iterable(index_of_sv_control))

    index_grouped_sample = group_similar_indices(
        index_of_sv_sample, distance=5, min_coverage=max(5, coverage_sample * 0.05)
    )
    index_grouped_control = group_similar_indices(
        index_of_sv_control, distance=5, min_coverage=max(5, coverage_control * 0.05)
    )

    index_converter_sample = define_index_converter(index_grouped_sample)
    index_converter_control = define_index_converter(index_grouped_control)

    index_and_sv_size_sample = [
        get_index_and_sv_size(m["CSSPLIT"], sv_type, MUTATION_LOCI, index_converter_sample)
        for m in io.read_jsonl(path_midsv_sample)
    ]
    index_and_sv_size_control = [
        get_index_and_sv_size(m["CSSPLIT"], sv_type, MUTATION_LOCI, index_converter_control)
        for m in io.read_jsonl(path_midsv_control)
    ]

    all_sv_index = {
        key
        for record in index_and_sv_size_sample + index_and_sv_size_control
        if isinstance(record, dict)
        for key in record.keys()
    }

    # If there are no SVs, return None
    if all_sv_index == set():
        return

    X_sample = extract_features(index_and_sv_size_sample, all_sv_index)
    X_control = extract_features(index_and_sv_size_control, all_sv_index)

    #######################################################
    # Clustering
    #######################################################

    X = np.vstack([X_sample, X_control])
    labels = optimize_labels(X, coverage_sample, coverage_control)

    #######################################################
    # Call consensus
    #######################################################

    ###################################################
    # 各IndexにおけるSV sizeのコンセンサスを、コントロールのアレルに反映させる
    ###################################################

    sv_index_by_label = defaultdict(set)
    for label, record in zip(labels, index_and_sv_size_sample):
        sv_index_by_label[label] |= record.keys()

    index_and_sv_size_by_label = defaultdict(list)
    for label, record in zip(labels, index_and_sv_size_sample):
        sv_index = sv_index_by_label[label]
        sv_size = {i: 0 for i in sv_index}
        sv_size |= record
        index_and_sv_size_by_label[label].append(sv_size)

    ###################################################
    # SV sizeのコンセンサスを取得
    ###################################################
    consensus_index_and_sv_size_by_label = defaultdict(dict)
    for label, index_and_sv_size in index_and_sv_size_by_label.items():
        frequency_counts = defaultdict(lambda: defaultdict(int))

        # Count the frequencies
        for record in index_and_sv_size:
            for key, value in record.items():
                frequency_counts[key][value] += 1

        # Extract the most frequent value for each key
        consensus_index_and_sv_size_by_label[label] = {
            key: max(values, key=values.get) for key, values in frequency_counts.items()
        }
        consensus_index_and_sv_size_by_label: dict[int, dict[int, int]] = dict(consensus_index_and_sv_size_by_label)

    if sv_type == "deletion" or sv_type == "inversion":
        ###################################################
        # ControlのアレルにSVを反映させる
        ###################################################

        cstag_by_label = {}
        fasta_by_label = {}

        def modify_tags(midsv_tags: list[str], start: int, end: int, modify_func: callable) -> list[str]:
            return [tag if not (start <= i <= end) else modify_func(tag) for i, tag in enumerate(midsv_tags)]

        for label, consensus_index_and_sv_size in consensus_index_and_sv_size_by_label.items():
            midsv_tags_control = ["=" + s for s in list(FASTA_ALLELES["control"])]
            for start, sv_size in consensus_index_and_sv_size.items():
                if sv_size == 0:  # Skip if the SV size is 0
                    continue
                end = start + sv_size
                if sv_type == "deletion":
                    midsv_tags_control = [
                        tag if not (start <= i <= end) else "-" + tag[1:] for i, tag in enumerate(midsv_tags_control)
                    ]
                elif sv_type == "inversion":
                    midsv_tags_control[start : end + 1] = revcomp_cssplits(midsv_tags_control[start : end + 1])

            cstag_by_label[label] = convert_cssplits_to_cstag(midsv_tags_control)
            fasta_by_label[label] = cstag.to_sequence(cstag_by_label[label])

    else:  # sv_type == "insertion"
        cssplits_iter: Iterator[list[list[str]]] = (m["CSSPLIT"].split(",") for m in io.read_jsonl(path_midsv_sample))

        cssplits_by_label = defaultdict(list)
        for label, cssplits in zip(labels, cssplits_iter):
            cssplits_by_label[label].append(cssplits)

        def get_consensus_insertion_by_label(
            cssplits_by_label: dict[int, list[str]], consensus_index_and_sv_size_by_label: dict[int, int]
        ) -> dict[int, dict[int, str]]:
            consensus_insertion_by_label = defaultdict(dict)
            for label, cssplits in cssplits_by_label.items():
                index_and_sv_size = consensus_index_and_sv_size_by_label[label]
                for index, size in index_and_sv_size.items():
                    insertions = []
                    if size == 0:
                        continue
                    for midsv_tag in cssplits:
                        if not midsv_tag[index].startswith("+"):
                            continue
                        if midsv_tag[index].count("|") == size:
                            insertions.append(midsv_tag[index])
                    consensus_insertion = Counter(insertions).most_common()[0][0]
                    consensus_insertion_by_label[label] |= {index: consensus_insertion}

            return dict(consensus_insertion_by_label)

        consensus_insertion_by_label = get_consensus_insertion_by_label(
            cssplits_by_label, consensus_index_and_sv_size_by_label
        )

        ###################################################
        # ControlのアレルにSV (Insertion) を反映させる
        ###################################################

        cstag_by_label = {}
        fasta_by_label = {}

        for label, consensus_insertion in consensus_insertion_by_label.items():
            midsv_tags_control = ["=" + s for s in list(FASTA_ALLELES["control"])]
            for index, insertion in consensus_insertion.items():
                midsv_tags_control[index] = insertion

            cstag_by_label[label] = convert_cssplits_to_cstag(midsv_tags_control)
            fasta_by_label[label] = cstag.to_sequence(cstag_by_label[label])

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
