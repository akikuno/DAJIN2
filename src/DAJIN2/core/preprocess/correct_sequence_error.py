from __future__ import annotations

import re
import json
import random
import tempfile
import subprocess
import numpy as np
from pathlib import Path
from itertools import islice
from typing import Generator
from collections import Counter, defaultdict


def _generate_3mer(midsv_sample: Generator[dict[str, str]]) -> Generator[str]:
    for cssplits in (cs["CSSPLIT"].split(",") for cs in midsv_sample):
        cssplits_kmer = []
        for idx_seq, cs in enumerate(cssplits):
            if idx_seq == 0:
                two_prior = cs
                cssplits_kmer.append(cs)
                continue
            if idx_seq == 1:
                one_prior = cs
                cssplits_kmer.append(cs)
                continue
            cssplits_kmer.append(two_prior + "," + one_prior + "," + cs)
            two_prior, one_prior = one_prior, cs
        yield cssplits_kmer


def _get_path_of_transposed_matrix(TEMPDIR, cssplits_sample: str, chunksize: int = 100, sep: str = "\t") -> str:
    # Write cssplits
    with tempfile.NamedTemporaryFile(dir=TEMPDIR) as fp:
        temp_cssplits = fp.name
    with open(temp_cssplits, "w") as out:
        for cssplits in cssplits_sample:
            out.write(sep.join(cssplits) + "\n")
    # Transpose by subprocess
    with tempfile.NamedTemporaryFile(dir=TEMPDIR) as fp:
        temp_output = fp.name
    with tempfile.NamedTemporaryFile(dir=TEMPDIR) as fp:
        temp_transposed_chunk = fp.name
    with tempfile.NamedTemporaryFile(dir=TEMPDIR) as fp:
        temp_pasted = fp.name
    with tempfile.NamedTemporaryFile(dir=TEMPDIR) as fp:
        output_path = fp.name
    with open(temp_cssplits, "r") as f:
        while True:
            next_n_lines = [line.strip().split(sep) for line in islice(f, chunksize)]
            if not next_n_lines:
                break
            with open(temp_transposed_chunk, "w") as out:
                for line in zip(*next_n_lines):
                    out.write(sep.join(line) + "\n")
            if not Path(temp_output).exists():
                subprocess.call(f"mv {temp_transposed_chunk} {temp_output}", shell=True)
            else:
                subprocess.call(f"paste -d '{sep}' {temp_output} {temp_transposed_chunk} > {temp_pasted}", shell=True)
                subprocess.call(f"cat {temp_pasted} > {temp_output}", shell=True)
        subprocess.call(f"cat {temp_output} > {output_path}", shell=True)
    return output_path


def _read_cssplits_as_tsv(input_path: str, sep: str = "\t"):
    with open(input_path, "r") as f:
        for line in f:
            yield line.strip().split(sep)


def _sampling_cssplits(transposed_midsv, mutation_loci) -> Generator:
    for idx_seq, cssplits_3mers in enumerate(transposed_midsv):
        sampling_2mer = defaultdict(list)
        if idx_seq == 0 or idx_seq == 1:
            yield defaultdict(Counter)
            continue
        for cssplits_3mer in cssplits_3mers:
            two_prior, one_prior, cs = cssplits_3mer.split(",")
            if cs[0] in mutation_loci[idx_seq] or cs.startswith("=") or cs == "N" or re.search(r"a|c|g|t|n", cs):
                previous_2mer = two_prior + "," + one_prior
                sampling_2mer[previous_2mer].append(cs)
            one_prior, two_prior = two_prior, cs
        yield {mer: Counter(samp) for mer, samp in sampling_2mer.items()}


def _correct_errors_cssplits(transposed_midsv, mutation_loci, sampling_cssplits):
    random.seed(1)
    for idx_seq, (cssplits_3mers, sampling_cssplit) in enumerate(zip(transposed_midsv, sampling_cssplits)):
        if idx_seq == 0 or idx_seq == 1:
            yield [cs for cs in cssplits_3mers]
            continue
        cssplits_corrected = []
        for cssplits_3mer in cssplits_3mers:
            two_prior, one_prior, cs = cssplits_3mer.split(",")
            if not (cs[0] in mutation_loci[idx_seq] or cs.startswith("=") or cs == "N" or re.search(r"a|c|g|t|n", cs)):
                previous_2mer = f"{two_prior},{one_prior}"
                try:
                    sampling = sampling_cssplit[previous_2mer]
                except KeyError:
                    cs = "N"
                else:
                    cs = random.choices(*zip(*sampling.items()))[0]
            cssplits_corrected.append(cs)
        yield cssplits_corrected


def _combine_with_qname_cssplits(midsv_sample, cssplits_corrected):
    for samp, cssplits in zip(midsv_sample, cssplits_corrected):
        yield {"QNAME": samp["QNAME"], "CSSPLIT": ",".join(cssplits)}


def read_json(filepath) -> Generator[dict[str, str]]:
    with open(filepath, "r") as f:
        for line in f:
            yield json.loads(line)


###############################################################################
# main
###############################################################################


def correct_sequence_error(TEMPDIR, FASTA_ALLELES, SAMPLE_NAME, CONTROL_NAME, MUTATION_LOCI_ALLELES) -> None:
    for allele in FASTA_ALLELES:
        mutation_loci = MUTATION_LOCI_ALLELES[allele]
        for name in [SAMPLE_NAME, CONTROL_NAME]:
            filepath = Path(TEMPDIR, "midsv", f"{name}_{allele}.json")
            midsv_3mer = _generate_3mer(read_json(filepath))
            path_transposed = _get_path_of_transposed_matrix(TEMPDIR, midsv_3mer, chunksize=1000)
            sampling = _sampling_cssplits(_read_cssplits_as_tsv(path_transposed), mutation_loci)
            cssplits_corrected_ = _correct_errors_cssplits(
                _read_cssplits_as_tsv(path_transposed), mutation_loci, sampling
            )
            path_transposed_corrected = _get_path_of_transposed_matrix(TEMPDIR, cssplits_corrected_, chunksize=1000)
            qname_cssplits = _combine_with_qname_cssplits(
                read_json(filepath), _read_cssplits_as_tsv(path_transposed_corrected)
            )
            out_filepath = Path(TEMPDIR, "midsv_corrected", f"{name}_{allele}.json")
            with open(out_filepath, "wt", encoding="utf-8") as f:
                for data in qname_cssplits:
                    f.write(json.dumps(data) + "\n")


# from __future__ import annotations

# import random
# import re
# from collections import Counter, defaultdict
# from typing import Generator
# import json
# from pathlib import Path


# def _sampling_cssplits(midsv_sample: dict[list[str, str]], mutation_loci: dict[str, set[int]]) -> list[dict[Counter]]:
#     sampling = [defaultdict(list) for _ in range(len(mutation_loci))]
#     for cssplits in (cs["CSSPLIT"] for cs in midsv_sample):
#         for idx_seq, cs in enumerate(cssplits.split(",")):
#             if idx_seq == len(mutation_loci) - 1:
#                 break
#             if idx_seq == 0:
#                 one_prior = cs
#                 continue
#             if idx_seq == 1:
#                 two_prior = cs
#                 continue
#             if cs[0] in mutation_loci[idx_seq] or cs.startswith("=") or cs == "N" or re.search(r"a|c|g|t|n", cs):
#                 prev_cs = one_prior + "," + two_prior
#                 sampling[idx_seq][prev_cs].append(cs)
#             one_prior, two_prior = two_prior, cs
#     sampling_cssplits = []
#     for samp in sampling:
#         samp_counter = dict()
#         for key, val in samp.items():
#             val = Counter(val)
#             samp_counter.update({key: val})
#         sampling_cssplits.append(samp_counter)
#     return sampling_cssplits


# def _correct_errors(
#     midsv_sample: Generator, mutation_loci: list(set(str)), sampling_cssplits: list(defaultdict(Counter))
# ) -> Generator[dict[str, str]]:
#     random.seed(1)
#     for samp in midsv_sample:
#         qname = samp["QNAME"]
#         cssplits = samp["CSSPLIT"].split(",")
#         for idx_seq, cs in enumerate(cssplits):
#             if idx_seq == len(mutation_loci) - 1:
#                 break
#             if idx_seq == 0:
#                 two_prior = cs
#                 continue
#             if idx_seq == 1:
#                 one_prior = cs
#                 continue
#            if not (cs[0] in mutation_loci[idx_seq] or cs.startswith("=") or cs == "N" or re.search(r"a|c|g|t|n", cs)):
#                 previous_2mer = f"{two_prior},{one_prior}"
#                 sampling = sampling_cssplits[idx_seq][previous_2mer]
#                 if sampling == Counter():
#                     cssplits[idx_seq] = "N"
#                     continue
#                 cssplits[idx_seq] = random.choices(*zip(*sampling.items()))[0]
#             two_prior, one_prior = one_prior, cs
#         cssplits_joined = ",".join(cssplits)
#         yield {"QNAME": qname, "CSSPLIT": cssplits_joined}


# def read_json(filepath) -> Generator[dict[str, str]]:
#     with open(filepath, "r") as f:
#         for line in f:
#             yield json.loads(line)


# def write_midsv(filepath, data) -> None:
#     with open(filepath, "wt", encoding="utf-8") as f:
#         for d in data:
#             f.write(json.dumps(d) + "\n")


# ###############################################################################
# # main
# ###############################################################################


# def correct_sequence_error(TEMPDIR, FASTA_ALLELES, SAMPLE_NAME, CONTROL_NAME, MUTATION_LOCI_ALLELES) -> None:
#     for allele in FASTA_ALLELES:
#         filepath_sample = Path(TEMPDIR, "midsv", f"{SAMPLE_NAME}_{allele}.json")
#         filepath_control = Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_{allele}.json")
#         # Extract mutation loci
#         mutation_loci = MUTATION_LOCI_ALLELES[allele]
#         # Correct sequence errors
#         sampling_sample = _sampling_cssplits(read_json(filepath_sample), mutation_loci)
#         sampling_control = _sampling_cssplits(read_json(filepath_control), mutation_loci)
#         midsv_corected_sample = _correct_errors(read_json(filepath_sample), mutation_loci, sampling_sample)
#         midsv_corected_control = _correct_errors(read_json(filepath_control), mutation_loci, sampling_control)
#         # Output corrected midsv
#         filepath_sample = Path(TEMPDIR, "midsv_corrected", f"{SAMPLE_NAME}_{allele}.json")
#         filepath_control = Path(TEMPDIR, "midsv_corrected", f"{CONTROL_NAME}_{allele}.json")
#         write_midsv(filepath_sample, midsv_corected_sample)
#         write_midsv(filepath_control, midsv_corected_control)
