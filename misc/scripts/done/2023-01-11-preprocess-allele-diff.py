from __future__ import annotations

"""
mappyを用いて異なるFASTA_ALLELEの特徴的変異を抽出する 例：点変異部位
"""

import midsv
from collections import defaultdict
from pathlib import Path
from itertools import permutations
from src.DAJIN2.core.preprocess import mappy_align

SAMPLE, CONTROL, ALLELE, NAME, GENOME, DEBUG, THREADS = (
    "examples/pm-tyr/barcode31.fq.gz",
    "examples/pm-tyr/barcode32.fq.gz",
    "examples/pm-tyr/design_tyr.fa",
    "test-pm-tyr",
    "mm10",
    True,
    14,
)


TEMPDIR = Path("DAJINResults", ".tempdir", NAME)


def extract_diff_loci(TEMPDIR) -> defaultdict[dict]:
    """
    Extract differencial loci between alleles
        - The purpose is to lower match_score between very similar alleles such as point mutation.
    """
    fasta_alleles = list(Path(TEMPDIR, "fasta").iterdir())
    mutation_alleles = defaultdict(dict)
    for comb in list(permutations(fasta_alleles, 2)):
        ref, query = comb
        ref_allele = ref.stem
        alignments = mappy_align.to_sam(ref, query, preset="splice")
        alignments = list(alignments)
        alignments = [a.split("\t") for a in alignments]
        alignments_midsv = midsv.transform(alignments, midsv=False, cssplit=True, qscore=False)[0]
        cssplits = alignments_midsv["CSSPLIT"].split(",")
        mutations = dict()
        for i, cs in enumerate(cssplits):
            if cs.startswith("="):
                continue
            mutations.update({i: cs})
        if len(mutations) < 10:
            mutation_alleles[ref_allele].update(mutations)
    return mutation_alleles
