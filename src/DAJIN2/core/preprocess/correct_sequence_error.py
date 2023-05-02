from __future__ import annotations

import random
import re
from collections import defaultdict, Counter
from typing import Generator


def _sampling_cssplits(midsv_sample: dict[list[str, str]], mutation_loci: dict[str, set[int]]) -> list[dict[Counter]]:
    sampling = [defaultdict(list) for _ in range(len(mutation_loci))]
    for cssplits in (cs["CSSPLIT"].split(",") for cs in midsv_sample):
        for idx_seq, cs in enumerate(cssplits):
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


def _correct_errors(midsv_sample, mutation_loci, sampling_cssplits) -> Generator[dict[str, str]]:
    random.seed(1)
    for samp in midsv_sample:
        qname = samp["QNAME"]
        cssplits = samp["CSSPLIT"].split(",")
        for idx_seq, cs in enumerate(cssplits):
            if idx_seq == len(mutation_loci) - 1:
                break
            if idx_seq == 0:
                one_prior = cs
                continue
            if idx_seq == 1:
                two_prior = cs
                continue
            if not (cs[0] in mutation_loci[idx_seq] or cs.startswith("=") or cs == "N" or re.search(r"a|c|g|t|n", cs)):
                prev_cs = f"{one_prior},{two_prior}"
                try:
                    sampling = sampling_cssplits[idx_seq][prev_cs]
                    cssplits[idx_seq] = random.choices(*zip(*sampling.items()))[0]
                except KeyError:
                    cssplits[idx_seq] = "N"
            one_prior, two_prior = two_prior, cs
        cssplits_joined = ",".join(cssplits)
        yield {"QNAME": qname, "CSSPLIT": cssplits_joined}


###############################################################################
# main
###############################################################################


def correct_sequence_error(midsv_sample_alleles, midsv_control_alleles, FASTA_ALLELES, MUTATION_LOCI_ALLELES) -> None:
    midsv_alleles_corrected = defaultdict(dict)
    for allele, sequence in FASTA_ALLELES.items():
        midsv_sample = midsv_sample_alleles[allele]
        midsv_control = midsv_control_alleles[allele]
        # Extract mutation loci
        mutation_loci = MUTATION_LOCI_ALLELES[allele]
        # # Correct sequence errors
        sampling_sample =_sampling_cssplits(midsv_sample, mutation_loci)
        sampling_control =_sampling_cssplits(midsv_control, mutation_loci)
        midsv_corected_sample = _correct_errors(midsv_sample, mutation_loci, sampling_sample)
        midsv_corected_control = _correct_errors(midsv_control, mutation_loci, sampling_control)
        midsv_alleles_corrected["sample"][allele] = midsv_corected_sample
        midsv_alleles_corrected["control"][allele] = midsv_corected_control
    return midsv_alleles_corrected
