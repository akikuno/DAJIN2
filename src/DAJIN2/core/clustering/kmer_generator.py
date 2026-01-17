from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

from DAJIN2.utils import fileio


def compress_insertion(midsv: list[str], index: int, compress_ins: bool) -> str:
    """Process insertion mutations based on the compress_ins flag.
    compress_ins: True -> +I{inserted sequence}
    compress_ins: False -> +{number of |}{inserted sequence}

    Note: Clustering is based on the number of nucleotides inserted, but this creates too fine a cluster, so the number of nucleotides is ignored.Ideally, clustering should be possible without ignoring the number of bases.
    """
    if compress_ins:
        for j in range(index - 1, index + 2):
            if midsv[j].startswith("+") and not midsv[j].startswith("+I"):
                midsv[j] = f"+I{midsv[j].split('|')[-1]}"
    else:
        for j in range(index - 1, index + 2):
            if midsv[j].startswith("+") and not midsv[j][1].isdigit():
                midsv[j] = f"+{midsv[j].count('|')}{midsv[j].split('|')[-1]}"
    return midsv


def annotate_cs_tag(midsv: str, mutation: set[str]) -> str:
    cs = midsv[0]  # +, -, *, or =
    if cs in mutation:
        return midsv
    else:
        return "@"


def generate_mutation_kmers(
    path_sample: Path | str, mutation_loci: list[set[str]], compress_ins: bool = True
) -> Iterator[list[str]]:
    midsv_sample = fileio.read_jsonl(path_sample)
    for midsv in (cs["MIDSV"].split(",") for cs in midsv_sample):
        mutation_kmers = ["@,@,@"]
        for i in range(1, len(midsv) - 1):
            if mutation_loci[i] == set():
                mutation_kmers.append("@,@,@")
                continue
            cs_current = midsv[i][0]  # +, -, *, or =
            if cs_current == "+":
                midsv = compress_insertion(midsv, i, compress_ins)
            if cs_current in mutation_loci[i]:
                cs_prev = annotate_cs_tag(midsv[i - 1], mutation_loci[i - 1])
                cs_next = annotate_cs_tag(midsv[i + 1], mutation_loci[i + 1])
                kmer = ",".join([cs_prev, midsv[i], cs_next])
                mutation_kmers.append(kmer)
            else:
                mutation_kmers.append("@,@,@")
        mutation_kmers.append("@,@,@")
        yield mutation_kmers
