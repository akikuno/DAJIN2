from __future__ import annotations

from collections import defaultdict
from itertools import groupby
from pathlib import Path
import pickle
import midsv
import numpy as np
import rapidfuzz
from sklearn.cluster import MeanShift, MiniBatchKMeans


def extract_insertions(path_sample: str, path_control: str, mutation_loci: dict) -> defaultdict[dict[dict]]:
    insertion_sample = defaultdict(int)
    insertion_control = defaultdict(int)
    for m in midsv.read_jsonl(path_sample):
        cssplits = m["CSSPLIT"].split(",")
        for idx in (i for i, m in enumerate(mutation_loci) if "+" in m):
            if cssplits[idx].startswith("+"):
                insertion_sample[frozenset((idx, cssplits[idx]))] += 1
    for m in midsv.read_jsonl(path_control):
        cssplits = m["CSSPLIT"].split(",")
        for idx in (i for i, m in enumerate(mutation_loci) if "+" in m):
            if cssplits[idx].startswith("+"):
                insertion_control[frozenset((idx, cssplits[idx]))] += 1
    insertions = defaultdict(dict)
    for key in insertion_sample.keys() - insertion_control.keys():
        score = insertion_sample[key]
        idx, seq = key
        if isinstance(seq, int):
            idx, seq = seq, idx
        insertions[idx].update({seq: score})
    return insertions


def group_consecutive_insertions(mutation_loci: dict[str, set[int]]) -> list[tuple[int]]:
    index = sorted(i for i, m in enumerate(mutation_loci) if "+" in m)
    index_grouped = []
    for _, group in groupby(enumerate(index), lambda i_x: i_x[0] - i_x[1]):
        items = [value for _, value in group]
        index_grouped.append(tuple(items))
    return index_grouped


def group_insertions(insertions, index_grouped) -> dict[dict[str, int]]:
    """
    Insertion bases in consecutive insertion base positions are grouped together in mutation_loci,
    as the insertion site may shift by a few bases.
    """
    insertions_grouped = defaultdict(dict)
    for idx in index_grouped:
        for i in idx:
            insertions_grouped[idx].update(insertions[i])
    return insertions_grouped


def merge_similar_insertions(insertions_grouped) -> dict[dict[frozenset[str], int]]:
    """
    Inserted bases with a high degree of similarity in terms of sequence and sequence length are regarded as the same.
    """
    insertions_merged = defaultdict(dict)
    for idx in insertions_grouped:
        insertion = insertions_grouped[idx]
        query = list(insertion.keys())[0]
        choices = list(insertion.keys())
        sequences_insertion = []
        scores_insertion = []
        for res in rapidfuzz.process.extract_iter(query, choices):
            sequences_insertion.append(res[0])
            scores_insertion.append(res[1])
        X = np.array(scores_insertion).reshape(-1, 1)
        X = (X - X.min()) / (X.max() - X.min())
        # count the base number of insertions
        counts = np.array([s.count("|") for s in sequences_insertion]).reshape(-1, 1)
        X = np.concatenate([X, counts], axis=1)
        clustering = MiniBatchKMeans(n_clusters=2, n_init="auto").fit(X)
        labels = clustering.labels_
        # count scores
        scores_similar = []
        scores_dissimilar = []
        label_first = labels[0]
        for label, val in zip(labels, insertion.values()):
            if label == label_first:
                scores_similar.append(val)
            else:
                scores_dissimilar.append(val)
        scores_similar = sum(scores_similar)
        scores_dissimilar = sum(scores_dissimilar)
        # combine sequences and scores
        seq_similar = set()
        seq_dissimilar = set()
        for label, seq in zip(labels, insertion.keys()):
            if label == label_first:
                seq_similar.add(seq)
            else:
                seq_dissimilar.add(seq)
        seq_similar = frozenset(seq_similar)
        seq_dissimilar = frozenset(seq_dissimilar)
        insertions_merged[idx].update({**{seq_similar: scores_similar}, **{seq_dissimilar: scores_dissimilar}})
    return insertions_merged


def extract_score_and_sequence(path_sample, insertions_merged) -> list[tuple[list[int], str]]:
    scores = []
    sequences = []
    for m in midsv.read_jsonl(path_sample):
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
    X = (X - X.min()) / (X.max() - X.min())
    clustering = MeanShift().fit(X)
    labels = clustering.labels_
    return labels


###########################################################
# consensus calling
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


def _call_percentage(cssplits: list[str]) -> list[dict[str, float]]:
    """call position weight matrix in defferent loci."""
    coverage = len(cssplits)
    cssplits_transposed = [list(cs) for cs in zip(*cssplits)]
    cons_percentage = []
    for cs_transposed in cssplits_transposed:
        count_cs = defaultdict(float)
        for cs in cs_transposed:
            count_cs[cs] += 1 / coverage * 100
        cons_percentage.append(count_cs)
    return cons_percentage


def _call_sequence(cons_percentage: list[dict[str, float]]) -> str:
    consensus_sequence = []
    for cons_per in cons_percentage:
        cons = max(cons_per, key=cons_per.get)
        if cons.startswith("="):
            cons = cons.replace("=", "")
        elif cons.startswith("-"):
            continue
        elif cons.startswith("*"):
            cons = cons[-1]
        elif cons.startswith("+"):
            cons_ins = cons.split("|")
            if cons_ins[-1].startswith("="):
                cons = cons.replace("=", "")
            elif cons_ins[-1].startswith("-"):
                cons = "".join(cons_ins[:-1])
            elif cons_ins[-1].startswith("*"):
                cons = "".join([*cons_ins[:-1], cons_ins[-1][-1]])
            cons = cons.replace("+", "")
            cons = cons.replace("|", "")
        consensus_sequence.append(cons)
    return ",".join(consensus_sequence)


def call_consensus(sequences_subset: list[dict]) -> dict[int, str]:
    cons_sequence = dict()
    sequences_subset.sort(key=lambda x: x["LABEL"])
    for label, group in groupby(sequences_subset, key=lambda x: x["LABEL"]):
        cssplits = [cs["CSSPLIT"].split(",") for cs in group]
        cons_per = _call_percentage(cssplits)
        cons_seq = _call_sequence(cons_per)
        cons_sequence[label] = cons_seq
    return cons_sequence


def remove_all_n(cons_sequence: dict[int, str]) -> dict[int, str]:
    cons_sequence_removed = dict()
    for label, seq in cons_sequence.items():
        if all(True if s == "N" else False for s in seq.split(",")):
            continue
        cons_sequence_removed[label] = seq
    return cons_sequence_removed


def extract_index_of_insertions(insertions, index_grouped) -> list[int]:
    index_of_insertions = []
    for idx_group in index_grouped:
        max_val = -1
        for idx in idx_group:
            if max_val < sum(insertions[idx].values()):
                max_val = sum(insertions[idx].values())
                max_idx = idx
        index_of_insertions.append(max_idx)
    return index_of_insertions


def generate_consensus(cons_sequence, index_of_insertions, sequence) -> dict[int, str]:
    consensus_sequence_insertion = dict()
    for label, cons_seq in cons_sequence.items():
        cons_seq = cons_seq.split(",")
        list_sequence = list(sequence)
        for idx, seq in zip(index_of_insertions, cons_seq):
            if seq == "N":
                continue
            list_sequence[idx] = seq
        consensus_sequence_insertion[label] = "".join(list_sequence)
    return consensus_sequence_insertion


###########################################################
# generate fasta
###########################################################


def update_labels(d: dict, prefix: str = "insertion") -> dict:
    keta = len(str(len(d)))
    return {f"{prefix}{i+1:0{keta}}": seq for i, seq in enumerate(d.values())}


def save_fasta(TEMPDIR: Path | str, consensus_sequence_insertion: dict) -> None:
    for header, seq in consensus_sequence_insertion.items():
        Path(TEMPDIR, "fasta", f"{header}.fasta").write_text(f">{header}\n{seq}")


###########################################################
# main
###########################################################


def generate_insertion_fasta(TEMPDIR, SAMPLE_NAME, CONTROL_NAME, FASTA_ALLELES) -> None:
    path_sample = Path(TEMPDIR, "midsv", f"{SAMPLE_NAME}_control.json")
    path_control = Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_control.json")
    with open(Path(TEMPDIR, "mutation_loci", f"{SAMPLE_NAME}_control.pickle"), "rb") as p:
        mutation_loci = pickle.load(p)
    sequence = FASTA_ALLELES["control"]
    insertions = extract_insertions(path_sample, path_control, mutation_loci)
    if not insertions:
        return None
    index_grouped = group_consecutive_insertions(mutation_loci)
    insertions_grouped = group_insertions(insertions, index_grouped)
    insertions_merged = merge_similar_insertions(insertions_grouped)
    insertions_scores_sequences = extract_score_and_sequence(path_sample, insertions_merged)
    labels = clustering_insertions([score for score, _ in insertions_scores_sequences])
    insertion_sequences_subset = subset_sequences([seq for _, seq in insertions_scores_sequences], labels, num=1000)
    cons_sequence = call_consensus(insertion_sequences_subset)
    cons_sequence_removed = remove_all_n(cons_sequence)
    index_of_insertions = extract_index_of_insertions(insertions, index_grouped)
    consensus_sequence_insertion = generate_consensus(cons_sequence_removed, index_of_insertions, sequence)
    consensus_sequence_insertion = update_labels(consensus_sequence_insertion)
    save_fasta(TEMPDIR, consensus_sequence_insertion)