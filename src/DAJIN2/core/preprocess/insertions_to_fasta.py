from __future__ import annotations

from pathlib import Path
from itertools import groupby
from collections import defaultdict
from typing import Generator

import numpy as np
from sklearn import metrics
from rapidfuzz import process
from rapidfuzz.distance import DamerauLevenshtein
from sklearn.cluster import MeanShift, MiniBatchKMeans

from DAJIN2.utils import io, config
from DAJIN2.utils.cssplits_handler import call_sequence

config.set_warnings_ignore()

import cstag

###########################################################
# Detect insertion sequences
###########################################################


def count_insertions(midsv_sample: Generator | str, mutation_loci: dict) -> dict[tuple[int, str], int]:
    # Extract insertion sequences that are more than 10 base pair length
    insertion_index_sequences = defaultdict(list)
    coverage = 0
    for m_sample in midsv_sample:
        coverage += 1
        cssplits = m_sample["CSSPLIT"].split(",")
        insertion_loci = (i for i, mut in enumerate(mutation_loci) if "+" in mut)
        for idx in insertion_loci:
            if cssplits[idx].startswith("+") and cssplits[idx].count("|") > 10:
                insertion_index_sequences[idx].append(cssplits[idx])

    insertion_counts = defaultdict(int)
    threshold = int(coverage * 0.5 / 100)
    for i, insertion_sequences in insertion_index_sequences.items():
        # Count similar insertion sequences
        query = insertion_sequences[0]
        seqs, scores, _ = zip(
            *process.extract_iter(query, insertion_sequences, scorer=DamerauLevenshtein.normalized_distance)
        )
        labels = [int(score * 10) for score in scores]
        insertion_labels_sequences = [(label, seq) for label, seq in zip(labels, insertion_sequences)]

        # Filter low frequency insertion sequences
        insertion_labels_sequences = sorted(insertion_labels_sequences, key=lambda x: x[0])
        for _, group in groupby(insertion_labels_sequences, key=lambda x: x[0]):
            group = list(group)
            if len(group) < threshold:
                continue
            seqs = (s for _, s in group)
            for seq in seqs:
                insertion_counts[(i, seq)] += 1

    return dict(insertion_counts)


def create_insertions_dict(
    insertion_sample: dict[tuple[int, str], int], insertion_control: dict[tuple[int, str], int]
) -> defaultdict[int, dict[str, int]]:
    """Create a dictionary of insertions that are present in the sample but not in the control."""
    insertions = defaultdict(dict)
    for key in insertion_sample.keys() - insertion_control.keys():
        score = insertion_sample[key]
        idx, seq = key
        if isinstance(seq, int):
            idx, seq = seq, idx
        insertions[idx][seq] = score
    return insertions


def extract_insertions(
    path_sample: Path | str, path_control: Path | str, mutation_loci: dict
) -> defaultdict[int, dict[str, int]]:
    insertion_sample = count_insertions(io.read_jsonl(path_sample), mutation_loci)
    insertion_control = count_insertions(io.read_jsonl(path_control), mutation_loci)
    return create_insertions_dict(insertion_sample, insertion_control)


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


def _get_normalized_scores(query: str, choices: list[str]) -> np.ndarray:
    seqs, scores, _ = zip(*process.extract_iter(query, choices, scorer=DamerauLevenshtein.normalized_distance))
    scores_np = np.array(scores).reshape(-1, 1)
    normalized_scores = (scores_np - scores_np.min()) / (scores_np.max() - scores_np.min())
    counts = np.array([s.count("|") for s in seqs]).reshape(-1, 1)
    X = np.concatenate([normalized_scores, counts], axis=1)
    return X


def _optimize_labels(X: np.array) -> list[int]:
    sample_size = X.shape[0]
    labels_prev = list(range(sample_size))
    for i in range(1, sample_size):
        np.random.seed(seed=1)
        labels_current = MiniBatchKMeans(n_clusters=i, random_state=1, n_init="auto").fit_predict(X).tolist()
        silhuette = metrics.silhouette_score(X, labels_current) if i > 1 else 0
        mutual_info = metrics.adjusted_mutual_info_score(labels_prev, labels_current)
        # print(i, Counter(labels_current), round(silhuette, 2), round(mutual_info, 2) ) # ! DEBUG
        if i == 2 and silhuette < 0:
            return [1] * len(labels_current)
        if 0.9 < silhuette or 0.9 < mutual_info:
            return labels_current
        labels_prev = labels_current
    return labels_current


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
        query = list(insertion.keys())[0]
        X = _get_normalized_scores(query, list(insertion.keys()))
        labels = _optimize_labels(X)
        insertions_merged[idx] = _get_merged_insertion(insertion, labels)
    return insertions_merged


###########################################################
# Cluster insertion alleles
###########################################################


def extract_score_and_sequence(path_sample, insertions_merged) -> list[tuple[list[int], str]]:
    scores = []
    sequences = []
    for m in io.read_jsonl(path_sample):
        score = defaultdict(int)
        seq = defaultdict(lambda: "N")
        cssplits = m["CSSPLIT"].split(",")
        for idx_grouped, ins in insertions_merged.items():
            score[idx_grouped] = 0
            seq[idx_grouped] = "N"
            for idx in idx_grouped:
                for seqs, value in ins.items():
                    if cssplits[idx] in seqs:
                        score[idx_grouped] = value
                        seq[idx_grouped] = cssplits[idx]
        if any(score.values()):
            scores.append([s for _, s in sorted(score.items())])
            sequences.append(",".join(s for _, s in sorted(seq.items())))
    return [(score, sequence) for score, sequence in zip(scores, sequences)]


def clustering_insertions(scores) -> list[int]:
    X = np.array(scores)
    if X.max() - X.min() == 0:
        return [1] * len(X)
    X = (X - X.min()) / (X.max() - X.min())
    np.random.seed(seed=1)
    clustering = MeanShift().fit(X)
    labels = clustering.labels_
    return labels.tolist()


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
    user_defined_alleles = set(FASTA_ALLELES)
    d_values = list(d.values())
    len_d = len(d_values)
    digits_up = 0
    # Update labels to avoid duplicating user-specified alleles
    # (insertion1 -> insertion01 -> insertion001...)
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


def extract_index_of_insertions(insertions, insertions_merged) -> list[int]:
    index_of_insertions = []
    for idx_group in insertions_merged:
        max_val = -1
        for idx in idx_group:
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
    index_insertions = set(i for i, m in enumerate(mutation_loci) if "+" in m)
    if insertions.keys() & index_insertions == set():
        return None
    insertions_merged = merge_similar_insertions(insertions, mutation_loci)
    insertions_scores_sequences = extract_score_and_sequence(path_sample, insertions_merged)
    labels = clustering_insertions([score for score, _ in insertions_scores_sequences])
    labels_filtered, insertion_scores_sequences_filtered = filter_minor_label(
        path_sample, labels, insertions_scores_sequences, threshold=0.5
    )
    insertion_sequences_subset = subset_sequences(
        [seq for _, seq in insertion_scores_sequences_filtered], labels_filtered, num=1000
    )
    consensus_of_insertions = call_consensus_of_insertion(insertion_sequences_subset)
    if consensus_of_insertions == dict():
        # If there is no insertion sequence, return None
        # It is possible when all insertion sequence annotated as `N` that is filtered out
        return None
    consensus_of_insertions = update_labels(consensus_of_insertions, FASTA_ALLELES)
    index_of_insertions = extract_index_of_insertions(insertions, insertions_merged)
    fasta_insertions = generate_fasta(consensus_of_insertions, index_of_insertions, sequence)
    cstag_insertions = generate_cstag(consensus_of_insertions, index_of_insertions, sequence)
    save_fasta(TEMPDIR, SAMPLE_NAME, fasta_insertions)
    save_html(TEMPDIR, SAMPLE_NAME, cstag_insertions)
