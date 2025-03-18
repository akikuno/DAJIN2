from __future__ import annotations

from collections import Counter, defaultdict
from itertools import chain
from pathlib import Path

import numpy as np

from DAJIN2.core.clustering.clustering import optimize_labels
from DAJIN2.core.preprocess.sv_handler import add_unique_allele_keys, extract_unique_sv, save_fasta, save_midsv
from DAJIN2.utils import io
from DAJIN2.utils.cssplits_handler import convert_cssplits_to_sequence

###############################################################################
# Extract features for clustering
# Features: The start index of the SV and the size of the SV
###############################################################################


def get_index_of_sv(misdv_string: str, sv_type: str, mutation_loci: list[set[str]]) -> list[int]:
    """Return a list with the start index of SV using the MIDSV tags"""

    tags = misdv_string.split(",")
    index_of_sv = []
    previous_tag = tags[0]

    for i, current_tag in enumerate(tags):
        if sv_type in ("deletion", "inversion"):
            if sv_type == "deletion":
                current_tag_is_sv = current_tag.startswith("-") and "-" in mutation_loci[i]
                previous_tag_is_not_sv = not previous_tag.startswith("-")
            else:  # inversion
                current_tag_is_sv = current_tag.islower()
                previous_tag_is_not_sv = not previous_tag.islower()

            # If the current tag is an SV and the previous tag is not an SV
            if current_tag_is_sv and previous_tag_is_not_sv:
                index_of_sv.append(i)

        elif sv_type == "insertion":
            if current_tag.startswith("+") and "+" in mutation_loci[i]:
                index_of_sv.append(i)

        previous_tag = current_tag
    return index_of_sv


def group_similar_indices(index_of_sv: list[int], distance: int = 5, min_coverage: float = 5.0) -> list[list[int]]:
    """The presence of indels causes shifts in the start index of SVs detected by Nanopore.
    To correct these shifts, group SV start indices that are within a specified distance of each other."""
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
    misdv_string: str, sv_type: str, mutation_loci: list[set[str]], index_converter: dict[int, int] = None
) -> dict[int, int]:
    """Return a dictionary with the start index and SV size"""

    tags = misdv_string.split(",")
    result = {}
    i = 0

    while i < len(tags):
        sv_size = 0
        start_index = None

        if sv_type == "deletion":
            while i < len(tags) and tags[i].startswith("-") and "-" in mutation_loci[i]:
                sv_size += 1
                i += 1

        elif sv_type == "inversion":
            while i < len(tags) and tags[i].islower():
                sv_size += 1
                i += 1

        elif sv_type == "insertion" and tags[i].startswith("+") and "+" in mutation_loci[i]:
            start_index = i
            sv_size = tags[i].count("|")

        if sv_size < 10:
            i += 1
            continue

        start_index = start_index if start_index is not None else i - sv_size

        if index_converter is None:
            result |= {start_index: sv_size}
        elif start_index in index_converter:
            start_index = index_converter[start_index]
            result |= {start_index: sv_size}

        i += 1

    return result


def extract_features(index_and_sv_size: list[dict[int, int]], all_sv_index: set[str]) -> np.ndarray:
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


# Data Loading
def load_data(tempdir, sample_name, control_name):
    mutation_loci = io.load_pickle(Path(tempdir, sample_name, "mutation_loci", "control", "mutation_loci.pickle"))
    path_midsv_sample = Path(tempdir, sample_name, "midsv", "control", f"{sample_name}.jsonl")
    path_midsv_control = Path(tempdir, control_name, "midsv", "control", f"{control_name}.jsonl")
    coverage_sample = io.count_newlines(path_midsv_sample)
    coverage_control = io.count_newlines(path_midsv_control)

    return mutation_loci, path_midsv_sample, path_midsv_control, coverage_sample, coverage_control


def process_sv_indices(path_midsv, sv_type, mutation_loci, coverage) -> dict[int, int]:
    """Process SV indices and group them based on proximity."""
    index_of_sv = [
        list(get_index_and_sv_size(m["CSSPLIT"], sv_type, mutation_loci, None).keys())
        for m in io.read_jsonl(path_midsv)
    ]
    index_of_sv = sorted(chain.from_iterable(index_of_sv))
    grouped_indices = group_similar_indices(index_of_sv, distance=5, min_coverage=max(5, coverage * 0.05))

    return define_index_converter(grouped_indices)


def extract_sv_features(midsv_path, sv_type, mutation_loci, index_converter) -> list[dict[int, int]]:
    """Extract SV start indices and sizes."""
    return [
        get_index_and_sv_size(m["CSSPLIT"], sv_type, mutation_loci, index_converter) for m in io.read_jsonl(midsv_path)
    ]


def get_sv_index_by_label(labels, index_and_sv_size_sample) -> dict[int, list[dict[int, int]]]:
    sv_index_by_label = defaultdict(set)
    for label, record in zip(labels, index_and_sv_size_sample):
        sv_index_by_label[label] |= record.keys()

    index_and_sv_size_by_label = defaultdict(list)
    for label, record in zip(labels, index_and_sv_size_sample):
        sv_index = sv_index_by_label[label]
        sv_size = {i: 0 for i in sv_index}
        sv_size |= record
        index_and_sv_size_by_label[label].append(sv_size)

    return dict(index_and_sv_size_by_label)


def remove_invalid_sv(midsv_by_label: dict[int, list[str]]) -> dict[int, list[str]]:
    """
    If the first or last index is deletion (starts with '-'), exclude the SV.
    """
    return {
        label: tag for label, tag in midsv_by_label.items() if not (tag[0].startswith("-") or tag[-1].startswith("-"))
    }


def get_midsv_consensus_by_label(
    index_and_sv_size_by_label, path_midsv_sample, labels, sv_type, fasta_control
) -> dict[int, list[str]]:
    # Retrieve the most frequent index and SV size
    consensus_index_and_sv_size_by_label = defaultdict(dict)
    for label, index_and_sv_size in index_and_sv_size_by_label.items():
        frequency_counts = defaultdict(lambda: defaultdict(int))

        # Count the frequencies of the SV sizes
        for record in index_and_sv_size:
            for index, sv_size in record.items():
                frequency_counts[index][sv_size] += 1

        consensus_index_and_sv_size_by_label[label] = {
            key: max(values, key=values.get) for key, values in frequency_counts.items()
        }

    midsv_by_label = {}

    if sv_type == "deletion" or sv_type == "inversion":
        ###################################################
        # Reflect SV (deletion / inversion) in the Control allele
        # ###################################################

        for label, consensus_index_and_sv_size in consensus_index_and_sv_size_by_label.items():
            midsv_tags_control = ["=" + s for s in list(fasta_control)]
            for start, sv_size in consensus_index_and_sv_size.items():
                if sv_size == 0:  # Skip if the SV size is 0
                    continue
                end = start + sv_size - 1
                if sv_type == "deletion":
                    midsv_tags_control = [
                        tag if not (start <= i <= end) else "-" + tag[1:] for i, tag in enumerate(midsv_tags_control)
                    ]
                elif sv_type == "inversion":
                    midsv_tags_control[start : end + 1] = [tag.lower() for tag in midsv_tags_control[start : end + 1]]

            midsv_by_label[label] = midsv_tags_control

    else:  # sv_type == "insertion"

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
                    if insertions:
                        consensus_insertion = Counter(insertions).most_common()[0][0]
                        consensus_insertion_by_label[label] |= {index: consensus_insertion}

            return dict(consensus_insertion_by_label)

        cssplits_iter = (m["CSSPLIT"].split(",") for m in io.read_jsonl(path_midsv_sample))
        cssplits_by_label = defaultdict(list)
        for label, cssplits in zip(labels, cssplits_iter):
            cssplits_by_label[label].append(cssplits)

        consensus_insertion_by_label = get_consensus_insertion_by_label(
            cssplits_by_label, consensus_index_and_sv_size_by_label
        )

        ###################################################
        # Reflect SV (Insertion) in the Control allele
        # ###################################################

        for label, consensus_insertion in consensus_insertion_by_label.items():
            midsv_tags_control = ["=" + s for s in list(fasta_control)]
            for index, insertion in consensus_insertion.items():
                midsv_tags_control[index] = insertion

            midsv_by_label[label] = midsv_tags_control

    return remove_invalid_sv(midsv_by_label)


###############################################################################
# main
###############################################################################


def detect_sv_alleles(TEMPDIR: Path, SAMPLE_NAME: str, CONTROL_NAME: str, FASTA_ALLELES: dict, sv_type: str) -> None:
    """Detect structural variant alleles and reflect consensus on control alleles."""

    mutation_loci, path_midsv_sample, path_midsv_control, coverage_sample, coverage_control = load_data(
        TEMPDIR, SAMPLE_NAME, CONTROL_NAME
    )

    #######################################################
    # Extract features: SV start index and SV size
    #######################################################

    index_converter_sample = process_sv_indices(path_midsv_sample, sv_type, mutation_loci, coverage_sample)
    if index_converter_sample == {}:
        return

    index_converter_control = process_sv_indices(path_midsv_control, sv_type, mutation_loci, coverage_control)

    index_and_sv_size_sample = extract_sv_features(path_midsv_sample, sv_type, mutation_loci, index_converter_sample)
    index_and_sv_size_control = extract_sv_features(
        path_midsv_control, sv_type, mutation_loci, index_converter_control
    )

    all_sv_indices = {
        key
        for record in index_and_sv_size_sample + index_and_sv_size_control
        if isinstance(record, dict)
        for key in record.keys()
    }

    # If there are no SVs, return None
    if all_sv_indices == set():
        return

    #######################################################
    # Clustering
    #######################################################

    X_sample = extract_features(index_and_sv_size_sample, all_sv_indices)
    X_control = extract_features(index_and_sv_size_control, all_sv_indices)
    X = np.vstack([X_sample, X_control])
    labels = optimize_labels(X, coverage_sample, coverage_control)

    threshold = max(coverage_sample * 0.05, 5)
    label_converter = {label: label if count > threshold else -1 for label, count in Counter(labels).items()}
    labels = [label_converter[label] for label in labels]

    #######################################################
    # Call consensus
    #######################################################

    # Aggregate the consensus of SV indices and sizes
    index_and_sv_size_by_label = get_sv_index_by_label(labels, index_and_sv_size_sample)

    # Consensus call of SVs
    midsv_by_label = get_midsv_consensus_by_label(
        index_and_sv_size_by_label, path_midsv_sample, labels, sv_type, FASTA_ALLELES["control"]
    )

    fasta_by_label = {label: convert_cssplits_to_sequence(midsv_tag) for label, midsv_tag in midsv_by_label.items()}

    # Discard -1 due to minor allele
    midsv_by_label = {label: tag for label, tag in midsv_by_label.items() if label != -1}
    fasta_by_label = {label: tag for label, tag in fasta_by_label.items() if label != -1}

    #######################################################
    # Remove similar alleles to user's alleles, or clustered alleles
    #######################################################

    fasta_by_label = extract_unique_sv(fasta_by_label, FASTA_ALLELES)
    midsv_by_label = {label: cs_tag for label, cs_tag in midsv_by_label.items() if label in fasta_by_label}

    #######################################################
    # Output cstag and fasta
    #######################################################
    midsv_by_label = add_unique_allele_keys(midsv_by_label, FASTA_ALLELES, key=sv_type)
    fasta_by_label = add_unique_allele_keys(fasta_by_label, FASTA_ALLELES, key=sv_type)

    save_midsv(TEMPDIR, SAMPLE_NAME, midsv_by_label)
    save_fasta(TEMPDIR, SAMPLE_NAME, fasta_by_label)
