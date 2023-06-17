from __future__ import annotations

from typing import Generator
import json
from pathlib import Path


def read_json(filepath: Path | str) -> Generator[dict[str, str]]:
    with open(filepath, "r") as f:
        for line in f:
            yield json.loads(line)


def generate_mutation_kmers(path_sample: Path | str, mutation_loci: list[set[str]]) -> Generator[list[str]]:
    midsv_sample = read_json(path_sample)
    for cssplit in (cs["CSSPLIT"].split(",") for cs in midsv_sample):
        mutation_kmers = ["N,N,N"]
        for i in range(1, len(cssplit) - 1):
            if mutation_loci[i] == set():
                mutation_kmers.append("N,N,N")
                continue
            mutation = cssplit[i][0]  # +, - , *, =, N
            if mutation == "+":
                cssplit[i] = "+I" + cssplit[i].split("|")[-1]
            if mutation in mutation_loci[i]:
                kmer = ",".join([cssplit[i - 1], cssplit[i], cssplit[i + 1]])
                mutation_kmers.append(kmer)
            else:
                mutation_kmers.append("N,N,N")
        mutation_kmers.append("N,N,N")
        yield mutation_kmers
