from __future__ import annotations

from pathlib import Path
from typing import Generator

from DAJIN2.utils import io


def compress_insertion(cssplit: list[str], index: int, compress_ins: bool) -> str:
    """Process insertion mutations based on the compress_ins flag.
    compress_ins: True -> +I{inserted sequence}
    compress_ins: False -> +{number of |}{inserted sequence}

    Note: Clustering is based on the number of nucleotides inserted, but this creates too fine a cluster, so the number of nucleotides is ignored.Ideally, clustering should be possible without ignoring the number of bases.
    """
    if compress_ins:
        for j in range(index - 1, index + 2):
            if cssplit[j].startswith("+") and not cssplit[j].startswith("+I"):
                cssplit[j] = f"+I{cssplit[j].split('|')[-1]}"
    else:
        for j in range(index - 1, index + 2):
            if cssplit[j].startswith("+") and not cssplit[j][1].isdigit():
                cssplit[j] = f"+{cssplit[j].count('|')}{cssplit[j].split('|')[-1]}"
    return cssplit


def annotate_cs_tag(cssplit: str, mutation: set[str]) -> str:
    cs = cssplit[0]  # +, - , *, =, or N
    if cs in mutation:
        return cssplit
    else:
        return "@"


def generate_mutation_kmers(
    path_sample: Path | str, mutation_loci: list[set[str]], compress_ins: bool = True
) -> Generator[list[str]]:
    midsv_sample = io.read_jsonl(path_sample)
    for cssplit in (cs["CSSPLIT"].split(",") for cs in midsv_sample):
        mutation_kmers = ["@,@,@"]
        for i in range(1, len(cssplit) - 1):
            if mutation_loci[i] == set():
                mutation_kmers.append("@,@,@")
                continue
            cs_current = cssplit[i][0]  # +, - , *, =, N
            if cs_current == "+":
                cssplit = compress_insertion(cssplit, i, compress_ins)
            if cs_current in mutation_loci[i]:
                cs_prev = annotate_cs_tag(cssplit[i - 1], mutation_loci[i - 1])
                cs_next = annotate_cs_tag(cssplit[i + 1], mutation_loci[i + 1])
                kmer = ",".join([cs_prev, cssplit[i], cs_next])
                mutation_kmers.append(kmer)
            else:
                mutation_kmers.append("@,@,@")
        mutation_kmers.append("@,@,@")
        yield mutation_kmers
