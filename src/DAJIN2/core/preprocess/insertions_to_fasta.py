from __future__ import annotations

import uuid
import random
from pathlib import Path
from itertools import groupby
from collections import defaultdict, Counter
from typing import Generator

import numpy as np
from rapidfuzz import process
from rapidfuzz.distance import DamerauLevenshtein
from sklearn.cluster import MeanShift

from DAJIN2.core.preprocess.mapping import to_sam
from DAJIN2.utils import io, config
from DAJIN2.utils.cssplits_handler import convert_cssplits_to_cstag

config.set_warnings_ignore()

import cstag


def remove_non_alphabets(s: str) -> str:
    """Convert a cssplits to a plain DNA sequence."""
    return "".join([char for char in s if char.isalpha()])


###########################################################
# Cluster insertion alleles
###########################################################


def clustering_insertions(insertions_cssplit: list[str]) -> list[int]:
    seq_all = [remove_non_alphabets(seq) for seq in insertions_cssplit]
    query = seq_all[0]
    _, distances, _ = zip(*process.extract_iter(query, seq_all, scorer=DamerauLevenshtein.normalized_distance))

    # Add random values from 0 to 1
    distances = list(distances)
    rng = np.random.default_rng(1)
    distances.extend(rng.random(max(100, len(seq_all))))

    return MeanShift().fit(np.array(distances).reshape(-1, 1)).labels_.tolist()[: len(seq_all)]


###########################################################
# Detect insertion sequences
###########################################################


def extract_all_insertions(midsv_sample: Generator, mutation_loci: list[set[str]]) -> dict[int, list[str]]:
    """To extract insertion sequences of **10 base pairs or more** at each index."""
    insertion_index_sequences_control = defaultdict(list)
    for m_sample in midsv_sample:
        cssplits = m_sample["CSSPLIT"].split(",")
        insertion_loci = (i for i, mut in enumerate(mutation_loci) if "+" in mut)
        for idx in insertion_loci:
            if cssplits[idx].startswith("+") and cssplits[idx].count("|") > 10:
                insertion_index_sequences_control[idx].append(cssplits[idx])

    return dict(insertion_index_sequences_control)


def extract_enriched_insertions(
    insertions_sample: dict, insertions_control: dict, coverage_sample: int
) -> dict[int, dict[str, int]]:
    enriched_insertions = dict()
    threshold_sample = max(5, int(coverage_sample * 0.5 / 100))
    for i in insertions_sample:
        ins_sample = insertions_sample[i]
        ins_control = insertions_control.get(i, [])

        labels_all = clustering_insertions(ins_sample + ins_control)
        labels_sample = labels_all[: len(ins_sample)]
        labels_control = labels_all[len(ins_sample) : len(ins_sample) + len(ins_control)]

        labels_count_sample = dict(Counter(labels_sample))
        labels_count_control = dict(Counter(labels_control))

        to_delete = set()
        # To remove labels containing a high proportion of controls (5% or more, or 5 or more reads).
        threshold_control = max(5, int(len(labels_control) * 0.05))
        for label, count_control in labels_count_control.items():
            if count_control > threshold_control:
                to_delete.add(label)

        # To delete labels with sample coverage below 0.5% or fewer than 5 reads.
        for label, count_sample in labels_count_sample.items():
            if count_sample < threshold_sample:
                to_delete.add(label)

        for label in to_delete:
            if label in labels_count_sample:
                del labels_count_sample[label]

        if labels_count_sample == dict():
            continue

        # Count the remaining insertion sequences.
        enriched_insertions[i] = dict(
            Counter([ins for ins, label in zip(ins_sample, labels_sample) if label in labels_count_sample])
        )
    return enriched_insertions


def extract_insertions(
    path_sample: Path, path_control: Path, mutation_loci: list[set[str]]
) -> dict[int, dict[str, int]]:
    insertions_sample = extract_all_insertions(io.read_jsonl(path_sample), mutation_loci)
    insertions_control = extract_all_insertions(io.read_jsonl(path_control), mutation_loci)
    coverage_sample = io.count_newlines(path_sample)

    return extract_enriched_insertions(insertions_sample, insertions_control, coverage_sample)


###########################################################
# Merge similar insertion sequences
###########################################################


def group_index_by_consecutive_insertions(mutation_loci: list[set[str]]) -> list[tuple[int]]:
    index = sorted(i for i, m in enumerate(mutation_loci) if "+" in m)
    index_grouped = []
    for _, group in groupby(enumerate(index), lambda i_x: i_x[0] - i_x[1]):
        items = [value for _, value in group]
        index_grouped.append(tuple(items))
    return index_grouped


def group_insertions(
    insertions: dict[int, dict[str, int]], index_grouped: list[tuple[int]]
) -> dict[tuple[int], dict[str, int]]:
    """
    Insertion bases in consecutive insertion base positions are grouped together in mutation_loci,
    as the insertion site may shift by a few bases.
    """
    insertions_grouped = defaultdict(dict)
    for idx in index_grouped:
        for i in idx:
            if i not in insertions:
                continue
            insertions_grouped[idx].update(insertions[i])
    return dict(insertions_grouped)


def get_merged_insertion(insertion: dict[str, int], labels: np.ndarray) -> dict[tuple[str], int]:
    insertion_label = [(label, {seq: count}) for label, (seq, count) in zip(labels, insertion.items())]
    insertion_label.sort(key=lambda x: x[0])

    insertion_merged = dict()
    for _, group in groupby(insertion_label, key=lambda x: x[0]):
        group = [g[1] for g in group]
        sequences, counts = set(), 0
        for seq_count in group:
            seq, count = next(iter(seq_count.items()))
            sequences.add(seq)
            counts += count
        insertion_merged[tuple(sorted(sequences))] = counts

    return insertion_merged


def merge_similar_insertions(
    insertions: dict[int, dict[str, int]], mutation_loci: list[set[str]]
) -> dict[tuple(int), dict[tuple[str], int]]:
    index_grouped = group_index_by_consecutive_insertions(mutation_loci)
    insertions_grouped = group_insertions(insertions, index_grouped)
    insertions_merged = dict()
    for idx, insertion in insertions_grouped.items():
        if len(insertion) == 1:
            seq, count = next(iter(insertion.items()))
            insertions_merged[idx] = {tuple([seq]): count}
            continue
        labels = clustering_insertions(insertion)
        insertions_merged[idx] = get_merged_insertion(insertion, labels)

    return insertions_merged


###########################################################
# Pre- and post-process of clustering insertion alleles
###########################################################


def flatten_keys_to_set(dict_keys: list[tuple[int]]) -> set[int]:
    flattened_nums = set()
    for key_tuple in dict_keys:
        for num in key_tuple:
            flattened_nums.add(num)
    return flattened_nums


def subset_insertions(
    insertions_merged: dict[tuple(int), dict[tuple[str], int]], num_subset: int = 100
) -> dict[tuple(int), dict[tuple[str], int]]:
    """
    If the number of seqs exceeds `num_subset`, limit it to only the `num_subset` elements with fiexd random sampling.
    """
    random.seed(1)
    insertions_merged_subset = defaultdict(dict)
    for idx_grouped, insertions in insertions_merged.items():
        for cs_insertion, counts in insertions.items():
            if len(cs_insertion) > num_subset:
                sequences = random.sample(cs_insertion, num_subset)
            else:
                sequences = cs_insertion
            sequences = tuple(remove_non_alphabets(seq) for seq in sequences)
            insertions_merged_subset[idx_grouped][sequences] = counts

    return dict(insertions_merged_subset)


def extract_score_and_sequence(
    path_sample, insertions_merged: dict[tuple(int), dict[tuple[str], int]]
) -> list[tuple[list[int], str]]:
    scores = []
    sequences = []

    insertions_merged_subset = subset_insertions(insertions_merged, num_subset=100)
    set_keys = flatten_keys_to_set(insertions_merged.keys())

    cache_score = defaultdict(int)
    for m in io.read_jsonl(path_sample):
        cssplits = m["CSSPLIT"].split(",")
        if not any(True for i in set_keys if cssplits[i].startswith("+")):
            continue

        score = [0] * len(insertions_merged)
        sequence = ["N"] * len(insertions_merged)
        for i, (idx_grouped, insertions) in enumerate(insertions_merged.items()):
            for idx in idx_grouped:
                if not cssplits[idx].startswith("+"):
                    continue

                for seqs, count in insertions.items():
                    if cssplits[idx] in seqs:
                        score[i] = count
                        sequence[i] = cssplits[idx]
                        continue

                if (idx, cssplits[idx]) in cache_score:
                    score[i] = cache_score[(idx, cssplits[idx])]
                    sequence[i] = cssplits[idx]
                    continue

                flag_break = False
                for seqs, count in insertions_merged_subset[idx_grouped].items():
                    for _, distance, _ in process.extract_iter(
                        remove_non_alphabets(cssplits[idx]), seqs, scorer=DamerauLevenshtein.normalized_distance
                    ):
                        if distance < 0.1:
                            score[i] = count
                            sequence[i] = cssplits[idx]
                            flag_break = True
                            break
                    if flag_break:
                        cache_score[(idx, cssplits[idx])] = count
                        break

        if any(score):
            scores.append(score)
            sequences.append(",".join(sequence))

    return [(score, sequence) for score, sequence in zip(scores, sequences)]


def filter_minor_label(
    path_sample: str,
    labels: list[int],
    insertions_scores_sequences: list[tuple[list[int], str]],
    threshold: float = 0.5,
) -> tuple(list[int], list[str]):
    coverage = io.count_newlines(path_sample)
    _, counts = np.unique(labels, return_counts=True)
    minor_labels = {label for label, count in enumerate(counts) if count / coverage * 100 < threshold}
    index_minor_labels = {i for i, label in enumerate(labels) if label in minor_labels}
    labels_filtered = [label for i, label in enumerate(labels) if i not in index_minor_labels]
    score_seq_filterd = [
        score_seq for i, score_seq in enumerate(insertions_scores_sequences) if i not in index_minor_labels
    ]
    return labels_filtered, score_seq_filterd


###########################################################
# Call consensus
###########################################################


def subset_sequences(sequences, labels, num=1000) -> list[dict]:
    """
    Downsampling to the required number of sequences when there is a huge number of inserted sequences
    """
    sequences_subset = []
    tmp_sequences = []
    random.seed(1)
    for sequence, label in zip(sequences, labels):
        tmp_sequences.append({"CSSPLIT": sequence, "LABEL": label})
    tmp_sequences.sort(key=lambda x: x["LABEL"])
    for _, group in groupby(tmp_sequences, key=lambda x: x["LABEL"]):
        group = list(group)
        if len(group) > num:
            sequences_subset += random.sample(group, num)
        else:
            sequences_subset += group
    return sequences_subset


def get_cstag_position(sam_insertions: list[str]) -> tuple[list[str], list[int]]:
    cs_tags = []
    positions = []
    for sam in sam_insertions:
        if sam.startswith("@"):
            continue
        cs_tags.append(sam.split("\t")[-1])
        positions.append(int(sam.split("\t")[3]))
    return cs_tags, positions


def mapping_insertion(
    TEMPDIR: Path, SAMPLE_NAME: str, cs_transposed: list[str], consensus_length: int
) -> Generator[str]:
    # Temporarily cache the reference sequence
    cs_insertion = next(cs for cs in cs_transposed if cs != "N" and len(cs) == consensus_length)
    ref_seq = cstag.to_sequence(convert_cssplits_to_cstag([cs_insertion]))

    path_reference = Path(TEMPDIR, SAMPLE_NAME, "fasta", f"reference-{str(uuid.uuid4())}.fasta")
    Path(path_reference).write_text(f">reference_insertion\n{ref_seq}\n")

    # Temporarily cache the query sequences
    path_query = Path(TEMPDIR, SAMPLE_NAME, "fasta", f"query-{str(uuid.uuid4())}.fasta")
    with Path(path_query).open("w") as file:
        for i, cs_insertion in enumerate(cs_transposed):
            que_seq = cstag.to_sequence(convert_cssplits_to_cstag([cs_insertion]))
            file.write(f">query_insertion_{i}\n{que_seq}\n")

    # Switch the mappy settings to either long-read or short-read, depending on the length of the reference sequence
    if len(ref_seq) > 500:
        sam_insertions = to_sam(path_reference, path_query)
    else:
        options = {"k": 3, "min_dp_score": 1, "min_chain_score": 1, "min_cnt": 1}
        sam_insertions = to_sam(path_reference, path_query, preset="sr", options=options)

    return sam_insertions


def generate_consensus_insertions(TEMPDIR: Path, SAMPLE_NAME: str, cssplits: list[list[str]]) -> str:
    consensus_insertion = []
    cssplits_transposed = (list(cs) for cs in zip(*cssplits))
    for cs_transposed in cssplits_transposed:
        if all(True if cs == "N" else False for cs in cs_transposed):
            consensus_insertion.append("N")
            continue

        count_N = sum(1 for cs in cs_transposed if cs == "N")

        count_cs = defaultdict(int)
        for cs in cs_transposed:
            if cs == "N":
                continue
            count_cs[len(cs)] += 1

        consensus_length = max(count_cs, key=count_cs.get)
        count_consensus = sum(1 for cs in cs_transposed if len(cs) == consensus_length)

        if count_N > count_consensus:
            consensus_insertion.append("N")
            continue

        sam_insertions = mapping_insertion(TEMPDIR, SAMPLE_NAME, cs_transposed, consensus_length)

        cs_tags, positions = get_cstag_position(sam_insertions)
        if cs_tags == []:  # When not mapped at all, append `N`.
            consensus_insertion.append("N")
        else:
            cs_insertion_consensus = cstag.to_sequence(cstag.consensus(cs_tags, positions))
            cs_last = Counter(cs.split("|")[-1] for cs in cs_transposed if len(cs) == consensus_length).most_common()[
                0
            ][0]
            cssplits_insertion_consensus = "|".join("+" + cs for cs in cs_insertion_consensus[:-1]) + "|" + cs_last
            consensus_insertion.append(cssplits_insertion_consensus)

    return ",".join(consensus_insertion)


def remove_all_n(cons_sequence: dict[int, str]) -> dict[int, str]:
    cons_sequence_removed = dict()
    for label, seq in cons_sequence.items():
        if all(True if s == "N" else False for s in seq.split(",")):
            continue
        cons_sequence_removed[label] = seq
    return cons_sequence_removed


def update_labels(d: dict, FASTA_ALLELES: dict) -> dict:
    """
    Update labels to avoid duplicating user-specified alleles
    (insertion1 -> insertion01 -> insertion001...)
    """
    user_defined_alleles = set(FASTA_ALLELES)
    d_values = list(d.values())
    len_d = len(d_values)
    digits_up = 0
    while True:
        digits = len(str(len_d)) + digits_up
        d_updated = {f"insertion{i+1:0{digits}}": seq for i, seq in enumerate(d_values)}
        if user_defined_alleles.isdisjoint(set(d_updated)):
            break
        digits_up += 1
    return d_updated


def call_consensus_of_insertion(
    TEMPDIR: Path, SAMPLE_NAME: str, FASTA_ALLELES: dict, insertion_sequences_subset: list[dict]
) -> dict[int, str]:
    """Generate consensus cssplits."""
    consensus_insertion_cssplits = dict()
    insertion_sequences_subset.sort(key=lambda x: x["LABEL"])
    for label, group in groupby(insertion_sequences_subset, key=lambda x: x["LABEL"]):
        cssplits = [cs["CSSPLIT"].split(",") for cs in group]
        consensus_insertion_cssplits[label] = generate_consensus_insertions(TEMPDIR, SAMPLE_NAME, cssplits)

    return update_labels(remove_all_n(consensus_insertion_cssplits), FASTA_ALLELES)


###########################################################
# Generate cstag and FASTA
###########################################################


def extract_index_of_insertions(
    insertions: dict[int, dict[str, int]], insertions_merged: dict[tuple(int), dict[tuple[str], int]]
) -> list[int]:
    """`insertions_merged` contains multiple surrounding indices for a single insertion allele. Among them, select the one index where the insertion allele is most frequent."""
    index_of_insertions = []
    for idx_group in insertions_merged:
        max_val = -1
        for idx in idx_group:
            if idx not in insertions:
                continue
            if max_val < sum(insertions[idx].values()):
                max_val = sum(insertions[idx].values())
                max_idx = idx
        index_of_insertions.append(max_idx)
    return index_of_insertions


def generate_cstag(
    consensus_of_insertions: dict[str, str], index_of_insertions: list[int], sequence: str
) -> dict[str, str]:
    cstag_insertions = dict()
    for label, cons_seq in consensus_of_insertions.items():
        cons_seq = cons_seq.split(",")
        list_sequence = list(sequence)
        for idx, seq in zip(index_of_insertions, cons_seq):
            if seq == "N":
                continue
            list_sequence[idx] = convert_cssplits_to_cstag([seq]) + "="
        cstag_insertions[label] = "cs:Z:=" + "".join(list_sequence)
    return cstag_insertions


def generate_fasta(cstag_insertions: dict[str, str]) -> dict[str, str]:
    fasta_insertions = dict()
    for label, cs_tag in cstag_insertions.items():
        fasta_insertions[label] = cstag.to_sequence(cs_tag)
    return fasta_insertions


###########################################################
# Save cstag (HTML) and fasta
###########################################################


def save_fasta(TEMPDIR: Path | str, SAMPLE_NAME: str, fasta_insertions: dict[str, str]) -> None:
    for header, seq in fasta_insertions.items():
        Path(TEMPDIR, SAMPLE_NAME, "fasta", f"{header}.fasta").write_text(f">{header}\n{seq}\n")


def save_cstag(TEMPDIR: Path | str, SAMPLE_NAME: str, cstag_insertions: dict[str, str]) -> None:
    for header, cs_tag in cstag_insertions.items():
        Path(TEMPDIR, SAMPLE_NAME, "cstag", f"{header}.txt").write_text(cs_tag + "\n")


# def save_html(TEMPDIR: Path, SAMPLE_NAME: str, cstag_insertions: dict[int, str]) -> None:
#     for header, cs_tag in cstag_insertions.items():
#         html = cstag.to_html(cs_tag, f"{SAMPLE_NAME} {header}")
#         Path(TEMPDIR, "report", "HTML", SAMPLE_NAME, f"{header}.html").write_text(html)


###########################################################
# main
###########################################################


def generate_insertion_fasta(TEMPDIR, SAMPLE_NAME, CONTROL_NAME, FASTA_ALLELES) -> None:
    PATH_SAMPLE = Path(TEMPDIR, SAMPLE_NAME, "midsv", "control.json")
    PATH_CONTROL = Path(TEMPDIR, CONTROL_NAME, "midsv", "control.json")
    SEQUENCE = FASTA_ALLELES["control"]
    MUTATION_LOCI = io.load_pickle(Path(TEMPDIR, SAMPLE_NAME, "mutation_loci", "control.pickle"))

    insertions = extract_insertions(PATH_SAMPLE, PATH_CONTROL, MUTATION_LOCI)
    if insertions == dict():
        return None
    insertions_merged = merge_similar_insertions(insertions, MUTATION_LOCI)
    insertions_scores_sequences = extract_score_and_sequence(PATH_SAMPLE, insertions_merged)
    labels = clustering_insertions([cssplit for _, cssplit in insertions_scores_sequences])
    labels_filtered, insertion_scores_sequences_filtered = filter_minor_label(
        PATH_SAMPLE, labels, insertions_scores_sequences, threshold=0.5
    )
    insertion_sequences_subset = subset_sequences(
        [seq for _, seq in insertion_scores_sequences_filtered], labels_filtered, num=1000
    )
    consensus_of_insertions = call_consensus_of_insertion(
        TEMPDIR, SAMPLE_NAME, FASTA_ALLELES, insertion_sequences_subset
    )
    if consensus_of_insertions == dict():
        """
        If there is no insertion sequence, return None
        It is possible when all insertion sequence annotated as `N` that is filtered out
        """
        return None
    index_of_insertions = extract_index_of_insertions(insertions, insertions_merged)
    cstag_insertions = generate_cstag(consensus_of_insertions, index_of_insertions, SEQUENCE)
    fasta_insertions = generate_fasta(cstag_insertions)

    save_fasta(TEMPDIR, SAMPLE_NAME, fasta_insertions)
    save_cstag(TEMPDIR, SAMPLE_NAME, cstag_insertions)
    # save_html(TEMPDIR, SAMPLE_NAME, cstag_insertions)

    _ = [path.unlink() for path in Path(TEMPDIR, SAMPLE_NAME, "fasta").glob("reference-*.fasta")]
    _ = [path.unlink() for path in Path(TEMPDIR, SAMPLE_NAME, "fasta").glob("query-*.fasta")]
