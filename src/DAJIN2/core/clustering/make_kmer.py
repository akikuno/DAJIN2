from __future__ import annotations

from typing import Generator
import json
from pathlib import Path


def read_json(filepath: Path | str) -> Generator[dict[str, str]]:
    with open(filepath, "r") as f:
        for line in f:
            yield json.loads(line)


# def generate_mutation_kmers(path_sample: Path | str, mutation_loci: list[set], kmer: int = 5) -> Generator[list[str]]:
#     midsv_sample = read_json(path_sample)
#     kmer_n = ",".join(["N"] * kmer)
#     cssplits = (cs["CSSPLIT"].split(",") for cs in midsv_sample)
#     for cssplit in cssplits:
#         mutation_kmers = [kmer_n] * (kmer // 2)
#         for i in range((kmer // 2), len(cssplit) - (kmer // 2)):
#             if mutation_loci[i] == set():
#                 mutation_kmers.append(kmer_n)
#                 continue
#             mutations = [cssplit[j][0] for j in range(i - kmer + 1, i + 1)]
#             if "+" in mutations:
#                 for j in range(i - kmer + 1, i + 1):
#                     if cssplit[j][0] == "+":
#                         cssplit[j] = "+I" + cssplit[j].split("|")[-1]
#             if all(mutation in mutation_loci[j] for j, mutation in zip(range(i - kmer + 1, i + 1), mutations)):
#                 kmer = ",".join(cssplit[i - kmer + 1 : i + 1])
#                 mutation_kmers.append(kmer)
#             else:
#                 mutation_kmers.append(kmer_n)
#         mutation_kmers += [kmer_n] * (kmer // 2)
#         yield mutation_kmers


def generate_mutation_kmers(
    path_sample: Path | str, mutation_loci: list[set[str]], compress_ins: bool = True
) -> Generator[list[str]]:
    midsv_sample = read_json(path_sample)
    for cssplit in (cs["CSSPLIT"].split(",") for cs in midsv_sample):
        mutation_kmers = ["N,N,N"]
        for i in range(1, len(cssplit) - 1):
            if mutation_loci[i] == set():
                mutation_kmers.append("N,N,N")
                continue
            mutation = cssplit[i][0]  # +, - , *, =, N
            if mutation == "+":
                if compress_ins:
                    for j in range(i - 1, i + 2):
                        if cssplit[j].startswith("+") and not cssplit[j].startswith("+I"):
                            cssplit[j] = f"+I{cssplit[j].split('|')[-1]}"
                else:
                    for j in range(i - 1, i + 2):
                        if cssplit[j].startswith("+") and not cssplit[j][1].isdigit():
                            cssplit[j] = f"+{cssplit[j].count('|')}{cssplit[j].split('|')[-1]}"
            if mutation in mutation_loci[i]:
                kmer = ",".join([cssplit[i - 1], cssplit[i], cssplit[i + 1]])
                mutation_kmers.append(kmer)
            else:
                mutation_kmers.append("N,N,N")
        mutation_kmers.append("N,N,N")
        yield mutation_kmers
