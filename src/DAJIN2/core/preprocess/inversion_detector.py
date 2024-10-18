from __future__ import annotations

from collections.abc import Iterator
from itertools import groupby
from pathlib import Path

import numpy as np
from rapidfuzz import process
from rapidfuzz.distance import DamerauLevenshtein
from sklearn.cluster import MeanShift

from DAJIN2.core.consensus.consensus import call_percentage
from DAJIN2.core.preprocess.sv_handler import add_unique_allele_keys, extract_unique_sv, save_cstag, save_fasta
from DAJIN2.utils import config, io
from DAJIN2.utils.cssplits_handler import convert_cssplit_to_dna, convert_cssplits_to_cstag, revcomp_cssplits

config.set_warnings_ignore()

import cstag

###########################################################
# Detect insertion sequences
###########################################################


def extract_inversions(midsv_sample: Iterator[list[str]]) -> list[dict[int, int, str, str]]:
    """Extract sequence of inversions and their start and end indices."""
    inversions = []
    for m_sample in midsv_sample:
        cssplits = m_sample["CSSPLIT"].split(",")
        min_idx = float("inf")
        max_idx = -1
        cs_inversion = []
        for i, cs in enumerate(cssplits):
            if cs[-1].islower():
                min_idx = min(min_idx, i)
                max_idx = max(max_idx, i)
                cs_inversion.append(cs)
        if cs_inversion:
            inversions.append(
                {
                    "start": min_idx,
                    "end": max_idx,
                    "seq": convert_cssplit_to_dna(",".join(cs_inversion)),
                    "CSSPLIT": m_sample["CSSPLIT"],
                }
            )

    return inversions


def _convert_sequences_to_distances(sequences: list[str]) -> list[float]:
    """Convert sequences to distances using Damerau-Levenshtein distance."""
    query = sequences[0]
    _, distances, _ = zip(*process.extract_iter(query, sequences, scorer=DamerauLevenshtein.normalized_distance))
    return distances


def _convert_to_scores(start_end: Iterator[int, int], distances: list[float]) -> list[int, int, float]:
    scores = []
    for (start, end), dist in zip(start_end, distances):
        scores.append([start, end, dist])
    return scores


def clustering_inversions(inversions: list[dict[int, int, str, str]], coverage: int) -> list[int]:
    distances = _convert_sequences_to_distances([d["seq"] for d in inversions])
    scores = _convert_to_scores(((d["start"], d["end"]) for d in inversions), distances)

    X = np.array(scores)
    min_size = max(5, int(coverage * 0.5 / 100))

    return MeanShift(min_bin_freq=min_size, bin_seeding=True).fit_predict(X).tolist()


###########################################################
# Call consensus
###########################################################


def call_consensus_of_inversion(
    inversions: list[dict[int, int, str, str]], labels: list[int], mutation_loci: list[set[str]], sequence: str
) -> dict[int, dict[str, str]]:
    """Generate consensus cssplits."""

    consensus_of_inversions = {}

    # Append labels
    for i, d in enumerate(inversions):
        d["label"] = labels[i]

    inversions.sort(key=lambda x: x["label"])
    for label, group in groupby(inversions, key=lambda x: x["label"]):
        group = list(group)[:1000]  # Downsample to 1000 sequences
        cssplits = [cs["CSSPLIT"].split(",") for cs in group]
        cons_percentage = call_percentage(cssplits, mutation_loci, sequence)
        cons_cssplits = [max(cons_per, key=cons_per.get) for cons_per in cons_percentage]

        # Detect inversion regions
        min_idx = float("inf")
        max_idx = -1
        for i, cs in enumerate(cons_cssplits):
            if cs[-1].islower():
                min_idx = min(min_idx, i)
                max_idx = max(max_idx, i)

        cons_inversion = revcomp_cssplits(cons_cssplits[min_idx : max_idx + 1])
        cons_cssplits_inversion_corrected = cons_cssplits[:min_idx] + cons_inversion + cons_cssplits[max_idx + 1 :]

        cons_cstag = convert_cssplits_to_cstag(cons_cssplits_inversion_corrected)
        cons_sequence = []
        cons_sequence = cstag.to_sequence(cons_cstag)
        cons_cssplits_inversion_corrected = ",".join(cons_cssplits_inversion_corrected)

        consensus_of_inversions[label] = {
            "CSTAG": cons_cstag,
            "SEQ": cons_sequence,
            "CSSPLIT": cons_cssplits_inversion_corrected,
        }

    return consensus_of_inversions


###########################################################
# main
###########################################################


def detect_inversions(TEMPDIR, SAMPLE_NAME, CONTROL_NAME, FASTA_ALLELES) -> None:
    path_sample = Path(TEMPDIR, SAMPLE_NAME, "midsv", "control", f"{SAMPLE_NAME}.jsonl")
    # path_control = Path(TEMPDIR, CONTROL_NAME, "midsv", "control", f"{CONTROL_NAME}.jsonl")
    mutation_loci = io.load_pickle(Path(TEMPDIR, SAMPLE_NAME, "mutation_loci", "control", "mutation_loci.pickle"))
    sequence = FASTA_ALLELES["control"]

    inversions = extract_inversions(io.read_jsonl(path_sample))
    if not inversions:  # If there is no insertion, return None
        return None

    coverage: int = io.count_newlines(path_sample)

    labels = clustering_inversions(inversions, coverage)

    consensus_of_inversions = call_consensus_of_inversion(inversions, labels, mutation_loci, sequence)

    seq_inversions = {key: d["SEQ"] for key, d in consensus_of_inversions.items()}
    consensus_of_unique_inversions = {
        key: consensus_of_inversions[key] for key in extract_unique_sv(seq_inversions, FASTA_ALLELES)
    }
    if consensus_of_unique_inversions == {}:
        return None

    consensus_of_unique_inversions = add_unique_allele_keys(consensus_of_unique_inversions, FASTA_ALLELES, "inversion")

    save_fasta(TEMPDIR, SAMPLE_NAME, {key: d["SEQ"] for key, d in consensus_of_unique_inversions.items()})
    save_cstag(TEMPDIR, SAMPLE_NAME, {key: d["CSTAG"] for key, d in consensus_of_unique_inversions.items()})

    # ToDo: Use CSSPLIT to reflect inversion in HTML
    # ToDo: Subtract from Sample when Control also has inversion
