from __future__ import annotations

from itertools import permutations
from pathlib import Path

import midsv

from DAJIN2.core.preprocess import mappy_align


def extract_knockin_loci(TEMPDIR: str | Path) -> dict(set(int)):
    """
    Returns:
        dict(set): loci of knockin in each fasta pairs
    """
    fasta_alleles = list(Path(TEMPDIR, "fasta").iterdir())
    fasta_alleles = [f for f in fasta_alleles if f.suffix != ".fai"]
    KNOCKIN_LOCI_ALLELES = dict()
    for pair in permutations(fasta_alleles, 2):
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
        KNOCKIN_LOCI_ALLELES.update({ref_allele: knockin_loci})
    return KNOCKIN_LOCI_ALLELES
