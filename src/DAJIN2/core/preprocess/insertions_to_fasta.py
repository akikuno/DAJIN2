from __future__ import annotations

from pathlib import Path
from itertools import groupby
from collections import defaultdict, Counter
from typing import Generator

import numpy as np
from rapidfuzz import process
from rapidfuzz.distance import DamerauLevenshtein
from sklearn.cluster import MeanShift

from DAJIN2.utils import io, config
from DAJIN2.utils.cssplits_handler import call_sequence

config.set_warnings_ignore()

import cstag


def remove_non_alphabets(s):
    # Using list comprehension to keep only alphabet characters
    return "".join([char for char in s if char.isalpha()])


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


def extract_all_insertions(midsv_sample: Generator, mutation_loci: dict) -> dict[int, list[str]]:
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


def extract_insertions(path_sample: Path, path_control: Path, mutation_loci: dict) -> dict[int, dict[str, int]]:
    insertions_sample = extract_all_insertions(io.read_jsonl(path_sample), mutation_loci)
    insertions_control = extract_all_insertions(io.read_jsonl(path_control), mutation_loci)
    coverage_sample = io.count_newlines(path_sample)

    return extract_enriched_insertions(insertions_sample, insertions_control, coverage_sample)


###########################################################
# Merge similar insertion sequences
###########################################################


def _group_index_by_consecutive_insertions(mutation_loci: dict[str, set[int]]) -> list[tuple[int]]:
    index = sorted(i for i, m in enumerate(mutation_loci) if "+" in m)
    index_grouped = []
    for _, group in groupby(enumerate(index), lambda i_x: i_x[0] - i_x[1]):
        items = [value for _, value in group]
        index_grouped.append(tuple(items))
    return index_grouped


def _group_insertions(insertions, index_grouped) -> dict[dict[str, int]]:
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
    return insertions_grouped


def _get_merged_insertion(insertion: dict[str, int], labels: np.ndarray) -> dict[frozenset, int]:
    insertion_label = [(label, {key: val}) for label, (key, val) in zip(labels, insertion.items())]
    insertion_label.sort(key=lambda x: x[0])
    insertion_merged = dict()
    for _, group in groupby(insertion_label, key=lambda x: x[0]):
        group = [g[1] for g in group]
        seq, count = set(), 0
        for g in group:
            key, val = list(g.items())[0]
            seq.add(key)
            count += val
        insertion_merged[frozenset(seq)] = count
    return insertion_merged


def merge_similar_insertions(insertions, mutation_loci) -> dict[dict[frozenset[str], int]]:
    index_grouped = _group_index_by_consecutive_insertions(mutation_loci)
    insertions_grouped = _group_insertions(insertions, index_grouped)
    insertions_merged = defaultdict(dict)
    for idx, insertion in insertions_grouped.items():
        if len(insertion) == 1:
            key, val = list(insertion.items())[0]
            insertions_merged[idx][frozenset([key])] = val
            continue
        labels = clustering_insertions(insertion)
        insertions_merged[idx] = _get_merged_insertion(insertion, labels)
    return insertions_merged


###########################################################
# Cluster insertion alleles
###########################################################


def extract_score_and_sequence(path_sample, insertions_merged) -> list[tuple[list[int], str]]:
    scores = []
    sequences = []
    for m in io.read_jsonl(path_sample):
        cssplits = m["CSSPLIT"].split(",")
        score = [0] * len(insertions_merged)
        sequence = ["N"] * len(insertions_merged)
        for i, (idx_grouped, insertions) in enumerate(insertions_merged.items()):
            for idx in idx_grouped:
                if not cssplits[idx].startswith("+"):
                    continue
                for seqs, count in insertions.items():
                    _, distances, _ = zip(
                        *process.extract_iter(cssplits[idx], seqs, scorer=DamerauLevenshtein.normalized_distance)
                    )
                    if any(True for d in distances if d < 0.1):
                        score[i] = count
                        sequence[i] = cssplits[idx]
        if any(score):
            scores.append(score)
            sequences.append(",".join(sequence))
    return [(score, sequence) for score, sequence in zip(scores, sequences)]


def filter_minor_label(
    path_sample: str, labels: list[int], insertions_scores_sequences, threshold: float = 0.5
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


###########################################################
# Call consensus
###########################################################


def subset_sequences(sequences, labels, num=1000) -> list[dict]:
    """
    Downsampling to the required number of sequences when there is a huge number of inserted sequences
    """
    sequences_subset = []
    tmp_sequences = []
    for sequence, label in zip(sequences, labels):
        tmp_sequences.append({"CSSPLIT": sequence, "LABEL": label})
    tmp_sequences.sort(key=lambda x: x["LABEL"])
    for _, group in groupby(tmp_sequences, key=lambda x: x["LABEL"]):
        sequences_subset.extend(list(group)[:num])
    return sequences_subset


def _call_percentage(cssplits: list[list[str]]) -> list[dict[str, float]]:
    """call position weight matrix in defferent loci."""
    coverage = len(cssplits)
    cssplits_transposed = (list(cs) for cs in zip(*cssplits))
    cons_percentage = []
    for cs_transposed in cssplits_transposed:
        count_cs = defaultdict(float)
        for cs in cs_transposed:
            count_cs[cs] += 1 / coverage * 100
        cons_percentage.append(dict(count_cs))
    return cons_percentage


def _remove_all_n(cons_sequence: dict[int, str]) -> dict[int, str]:
    cons_sequence_removed = dict()
    for label, seq in cons_sequence.items():
        if all(True if s == "N" else False for s in seq.split(",")):
            continue
        cons_sequence_removed[label] = seq
    return cons_sequence_removed


def call_consensus_of_insertion(insertion_sequences_subset: list[dict]) -> dict[int, str]:
    cons_sequence = dict()
    insertion_sequences_subset.sort(key=lambda x: x["LABEL"])
    for label, group in groupby(insertion_sequences_subset, key=lambda x: x["LABEL"]):
        cssplits = [cs["CSSPLIT"].split(",") for cs in group]
        cons_per = _call_percentage(cssplits)
        cons_seq = call_sequence(cons_per, sep=",")
        cons_sequence[label] = cons_seq
    return _remove_all_n(cons_sequence)


def extract_index_of_insertions(
    insertions: dict[int, dict[str, int]], insertions_merged: dict[dict[frozenset[str], int]]
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


###########################################################
# generate and save fasta
###########################################################


def generate_fasta(cons_sequence, index_of_insertions, sequence) -> dict[int, str]:
    fasta_insertions = dict()
    for label, cons_seq in cons_sequence.items():
        cons_seq = cons_seq.split(",")
        list_sequence = list(sequence)
        for idx, seq in zip(index_of_insertions, cons_seq):
            if seq == "N":
                continue
            list_sequence[idx] = seq
        fasta_insertions[label] = "".join(list_sequence)
    return fasta_insertions


def save_fasta(TEMPDIR: Path | str, SAMPLE_NAME: str, fasta_insertions: dict) -> None:
    for header, seq in fasta_insertions.items():
        Path(TEMPDIR, SAMPLE_NAME, "fasta", f"{header}.fasta").write_text(f">{header}\n{seq}")


###########################################################
# generate and save as HTML and PDF
###########################################################


def generate_cstag(cons_sequence, index_of_insertions, sequence) -> dict[int, str]:
    cstag_insertions = dict()
    for label, cons_seq in cons_sequence.items():
        cons_seq = cons_seq.split(",")
        list_sequence = list(sequence)
        for idx, seq in zip(index_of_insertions, cons_seq):
            if seq == "N":
                continue
            seq = seq.lower()
            list_sequence[idx] = f"+{seq}="
        cstag_insertions[label] = "cs:Z:=" + "".join(list_sequence)
    return cstag_insertions


def save_html(TEMPDIR: Path, SAMPLE_NAME: str, cstag_insertions: dict) -> None:
    for header, cs_tag in cstag_insertions.items():
        html = cstag.to_html(cs_tag, f"{SAMPLE_NAME} {header}")
        Path(TEMPDIR, "report", "HTML", SAMPLE_NAME, f"{header}.html").write_text(html)


###########################################################
# main
###########################################################


def generate_insertion_fasta(TEMPDIR, SAMPLE_NAME, CONTROL_NAME, FASTA_ALLELES) -> None:
    path_sample = Path(TEMPDIR, SAMPLE_NAME, "midsv", "control.json")
    path_control = Path(TEMPDIR, CONTROL_NAME, "midsv", "control.json")
    sequence = FASTA_ALLELES["control"]
    mutation_loci = io.load_pickle(Path(TEMPDIR, SAMPLE_NAME, "mutation_loci", "control.pickle"))
    insertions = extract_insertions(path_sample, path_control, mutation_loci)
    if insertions == dict():
        return None
    insertions_merged = merge_similar_insertions(insertions, mutation_loci)
    insertions_scores_sequences = extract_score_and_sequence(path_sample, insertions_merged)
    labels = clustering_insertions([cssplit for _, cssplit in insertions_scores_sequences])
    labels_filtered, insertion_scores_sequences_filtered = filter_minor_label(
        path_sample, labels, insertions_scores_sequences, threshold=0.5
    )
    insertion_sequences_subset = subset_sequences(
        [seq for _, seq in insertion_scores_sequences_filtered], labels_filtered, num=1000
    )
    consensus_of_insertions = call_consensus_of_insertion(insertion_sequences_subset)
    if consensus_of_insertions == dict():
        """
        If there is no insertion sequence, return None
        It is possible when all insertion sequence annotated as `N` that is filtered out
        """
        return None
    consensus_of_insertions = update_labels(consensus_of_insertions, FASTA_ALLELES)
    index_of_insertions = extract_index_of_insertions(insertions, insertions_merged)
    fasta_insertions = generate_fasta(consensus_of_insertions, index_of_insertions, sequence)
    cstag_insertions = generate_cstag(consensus_of_insertions, index_of_insertions, sequence)
    save_fasta(TEMPDIR, SAMPLE_NAME, fasta_insertions)
    save_html(TEMPDIR, SAMPLE_NAME, cstag_insertions)
