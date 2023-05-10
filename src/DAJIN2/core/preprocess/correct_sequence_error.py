from __future__ import annotations

import random
import re
from collections import Counter, defaultdict
from typing import Generator
import json
from pathlib import Path


def _sampling_cssplits(midsv_sample: dict[list[str, str]], mutation_loci: dict[str, set[int]]) -> list[dict[Counter]]:
    sampling = [defaultdict(list) for _ in range(len(mutation_loci))]
    for cssplits in (cs["CSSPLIT"] for cs in midsv_sample):
        for idx_seq, cs in enumerate(cssplits.split(",")):
            if idx_seq == len(mutation_loci) - 1:
                break
            if idx_seq == 0:
                one_prior = cs
                continue
            if idx_seq == 1:
                two_prior = cs
                continue
            if cs[0] in mutation_loci[idx_seq] or cs.startswith("=") or cs == "N" or re.search(r"a|c|g|t|n", cs):
                prev_cs = one_prior + "," + two_prior
                sampling[idx_seq][prev_cs].append(cs)
            one_prior, two_prior = two_prior, cs
    sampling_cssplits = []
    for samp in sampling:
        samp_counter = dict()
        for key, val in samp.items():
            val = Counter(val)
            samp_counter.update({key: val})
        sampling_cssplits.append(samp_counter)
    return sampling_cssplits


def _correct_errors(
    midsv_sample: Generator, mutation_loci: list(set(str)), sampling_cssplits: list(defaultdict(Counter))
) -> Generator[dict[str, str]]:
    random.seed(1)
    for samp in midsv_sample:
        qname = samp["QNAME"]
        cssplits = samp["CSSPLIT"].split(",")
        for idx_seq, cs in enumerate(cssplits):
            if idx_seq == len(mutation_loci) - 1:
                break
            if idx_seq == 0:
                two_prior = cs
                continue
            if idx_seq == 1:
                one_prior = cs
                continue
            if not (cs[0] in mutation_loci[idx_seq] or cs.startswith("=") or cs == "N" or re.search(r"a|c|g|t|n", cs)):
                previous_2mer = f"{two_prior},{one_prior}"
                sampling = sampling_cssplits[idx_seq][previous_2mer]
                if sampling == Counter():
                    cssplits[idx_seq] = "N"
                    continue
                cssplits[idx_seq] = random.choices(*zip(*sampling.items()))[0]
            two_prior, one_prior = one_prior, cs
        cssplits_joined = ",".join(cssplits)
        yield {"QNAME": qname, "CSSPLIT": cssplits_joined}


def read_midsv(filepath) -> Generator[dict[str, str]]:
    with open(filepath, "r") as f:
        for line in f:
            yield json.loads(line)


def write_midsv(filepath, data) -> None:
    with open(filepath, "wt", encoding="utf-8") as f:
        for d in data:
            f.write(json.dumps(d) + "\n")


###############################################################################
# main
###############################################################################


def correct_sequence_error(TEMPDIR, FASTA_ALLELES, SAMPLE_NAME, CONTROL_NAME, MUTATION_LOCI_ALLELES) -> None:
    for allele in FASTA_ALLELES:
        filepath_sample = Path(TEMPDIR, "midsv", f"{SAMPLE_NAME}_{allele}.json")
        filepath_control = Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_{allele}.json")
        # Extract mutation loci
        mutation_loci = MUTATION_LOCI_ALLELES[allele]
        # Correct sequence errors
        sampling_sample = _sampling_cssplits(read_midsv(filepath_sample), mutation_loci)
        sampling_control = _sampling_cssplits(read_midsv(filepath_control), mutation_loci)
        midsv_corected_sample = _correct_errors(read_midsv(filepath_sample), mutation_loci, sampling_sample)
        midsv_corected_control = _correct_errors(read_midsv(filepath_control), mutation_loci, sampling_control)
        # Output corrected midsv
        filepath_sample = Path(TEMPDIR, "midsv_corrected", f"{SAMPLE_NAME}_{allele}.json")
        filepath_control = Path(TEMPDIR, "midsv_corrected", f"{CONTROL_NAME}_{allele}.json")
        write_midsv(filepath_sample, midsv_corected_sample)
        write_midsv(filepath_control, midsv_corected_control)
