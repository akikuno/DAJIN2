from __future__ import annotations

import re
from collections import Counter, defaultdict
from difflib import get_close_matches
from itertools import permutations
from pathlib import Path

import midsv
from scipy import stats
from scipy.spatial.distance import cosine

from DAJIN2.core.preprocess import mappy_align

###############################################################################
# functions
###############################################################################


def extract_knockin_loci(TEMPDIR) -> defaultdict(set):
    """
    Returns:
        defaultdict(set): loci of knockin in each fasta pairs
    """
    fasta_alleles = list(Path(TEMPDIR, "fasta").iterdir())
    fasta_alleles = [f for f in fasta_alleles if f.suffix != ".fai"]
    knockin_alleles = defaultdict(set)
    for pair in list(permutations(fasta_alleles, 2)):
        ref, query = pair
        ref_allele = ref.stem
        alignments = mappy_align.to_sam(ref, query, preset="splice")
        alignments = [a.split("\t") for a in alignments]
        alignments_midsv = midsv.transform(alignments, midsv=False, cssplit=True, qscore=False)[0]
        cssplits = alignments_midsv["CSSPLIT"].split(",")
        knockin_loci = set()
        for i, cs in enumerate(cssplits):
            if cs == "N" or cs.startswith("-"):
                knockin_loci.add(i)
        knockin_alleles[ref_allele] = knockin_loci
    return knockin_alleles


def get_5mer_of_knockin_loci(sequence: str, knockin_loci: list) -> dict:
    sequence_kmer = dict()
    for i in knockin_loci:
        sequence_kmer.update({i + 2: sequence[i : i + 5]})
    return sequence_kmer


def get_5mer_of_sequence(sequence: str) -> dict:
    if len(sequence) <= 5:
        return [sequence]
    sequence_kmer = dict()
    for i in range(len(sequence) - 5):
        sequence_kmer.update({i + 2: sequence[i : i + 5]})
    return sequence_kmer


def get_idx_of_similar_5mers(knockin_kmer: dict, sequence_kmer: dict, knockin_loci: set, n=100) -> defaultdict(set):
    idx_of_similar_5mers = defaultdict(list)
    for locus, kmer in knockin_kmer.items():
        idxes = set()
        for similar_kmer in get_close_matches(kmer, sequence_kmer.values(), n=n, cutoff=0.0):
            for idx in (key for key, val in sequence_kmer.items() if val == similar_kmer):
                if set(range(idx - 2, idx + 3)) & knockin_loci:
                    continue
                idxes.add(idx)
        idx_of_similar_5mers[locus] = idxes
    return idx_of_similar_5mers


def count_indel_5mer(cssplits_transposed: list, indexes: list(int)) -> defaultdict(dict):
    count_5mer = defaultdict(dict)
    for i in indexes:
        cssplits_5mer = cssplits_transposed[i - 2 : i + 3]
        count = {"ins": [1] * 5, "del": [1] * 5, "sub": [1] * 5}
        for j, cs in enumerate(cssplits_5mer):
            counter = Counter(cs)
            for key, cnt in counter.items():
                if key.startswith("=") or key == "N" or re.search(r"a|c|g|t|n", key):
                    continue
                if key.startswith("+"):
                    count["ins"][j] += cnt
                elif key.startswith("-"):
                    count["del"][j] += cnt
                elif key.startswith("*"):
                    count["sub"][j] += cnt
        count_5mer[i] = count
    return count_5mer


def replace_errors_to_match(cssplits_sample: list, sequence_errors: defaultdict(set), sequence: str):
    cssplits_replaced = []
    for cssplits in cssplits_sample:
        cssplits_copy = cssplits.copy()
        for i, error in sequence_errors.items():
            cssplits_5mer = cssplits_copy[i - 2 : i + 3]
            for j, mer in enumerate(cssplits_5mer):
                match_seq = "=" + sequence[i - 2 + j]
                if "ins" in error and mer.startswith("+"):
                    cssplits_5mer[j] = match_seq
                if "del" in error and mer.startswith("-"):
                    cssplits_5mer[j] = match_seq
                if "del" in error and mer.startswith("*"):
                    cssplits_5mer[j] = match_seq
            cssplits_copy[i - 2 : i + 3] = cssplits_5mer
        cssplits_replaced.append(cssplits_copy)
    return cssplits_replaced


##########################################################
# main
##########################################################


def execute(TEMPDIR, FASTA_ALLELES, CONTROL_NAME, SAMPLE_NAME):
    knockin_alleles = extract_knockin_loci(TEMPDIR)
    for allele, sequence in FASTA_ALLELES.items():
        sequence = FASTA_ALLELES[allele]
        knockin_loci = knockin_alleles[allele]
        midsv_sample = midsv.read_jsonl(Path(TEMPDIR, "midsv", f"{SAMPLE_NAME}_{allele}.jsonl"))
        midsv_control = midsv.read_jsonl(Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_{allele}.jsonl"))
        cssplits_sample = [m["CSSPLIT"].split(",") for m in midsv_sample]
        cssplits_control = [m["CSSPLIT"].split(",") for m in midsv_control]
        # Split the knock-in sequence and control into 5 mer and find 100 similar sequences in the control
        knockin_kmer = get_5mer_of_knockin_loci(sequence, knockin_loci)
        sequence_kmer = get_5mer_of_sequence(sequence)
        idx_of_similar_5mers = get_idx_of_similar_5mers(knockin_kmer, sequence_kmer, knockin_loci, n=100)
        # Find the number of indels in 5mer of similar sequence
        count_5mer_similar_sequences = defaultdict(dict)
        cssplits_transposed = [list(t) for t in zip(*cssplits_control)]
        for i, indexes in idx_of_similar_5mers.items():
            count_5mer_similar_sequences[i] = count_indel_5mer(cssplits_transposed, indexes)
        # Find the number of indels in 5mer of Knock-in sequence
        cssplits_transposed = [list(t) for t in zip(*cssplits_sample)]
        count_5mer_knockin = count_indel_5mer(cssplits_transposed, knockin_loci)
        # If there is an error profile similar to the Knock-in sequence and similar sequences, consider it a sequencing error
        coverage_sample = len(midsv_sample)
        coverage_control = len(midsv_control)
        sequence_errors = defaultdict(set)
        for i, count_knockin in count_5mer_knockin.items():
            knockin_mutation = defaultdict(list)
            for mutation in ["ins", "del", "sub"]:
                knockin_mutation[mutation] = [c / coverage_sample for c in count_knockin[mutation]]
            count_control = count_5mer_similar_sequences[i]
            for _, count in count_control.items():
                for mutation in ["ins", "del", "sub"]:
                    knockin = knockin_mutation[mutation]
                    control = [c / coverage_control for c in count[mutation]]
                    distance = 1 - cosine(knockin, control)
                    _, pvalue = stats.ttest_ind(knockin, control, equal_var=False)
                    if distance > 0.7 and pvalue > 0.01:
                        sequence_errors[i].add(mutation)
        # Correct sequencing errors in knock-in sequences
        cssplits_replaced = replace_errors_to_match(cssplits_sample, sequence_errors, sequence)
        # Replace CSSPLIT
        cssplits_corrected = [",".join(cs) for cs in cssplits_replaced]
        for i, cssplits in enumerate(cssplits_corrected):
            midsv_sample[i]["CSSPLIT"] = cssplits
        # Save to jsonl
        midsv.write_jsonl(midsv_sample, Path(TEMPDIR, "midsv", f"{SAMPLE_NAME}_{allele}.jsonl"))
