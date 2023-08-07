from __future__ import annotations

import midsv
import pickle
from pathlib import Path
from DAJIN2.core.preprocess import align


def select_control(path_fasta: list[Path]) -> Path:
    return [f for f in path_fasta if f.stem == "control"][0]


def select_others(path_fasta: list[Path]) -> list[Path]:
    return [f for f in path_fasta if f.stem != "control"]


def has_splice(alignments_midsv: dict[str]) -> bool:
    cssplits = alignments_midsv["CSSPLIT"].split(",")
    i = 0
    while i < len(cssplits):
        if cssplits[i].startswith("+") and (cssplits[i + 1].startswith("-") or cssplits[i + 1].startswith("N")):
            length_ins = cssplits[i].count("+")
            current_op = cssplits[i + 1][0]
            length_del = 0
            while current_op == cssplits[i + 1][0]:
                length_del += 1
                i += 1
            if length_ins == length_del:
                return True
        i += 1
    return False


def get_cssplit(reference: Path, query: Path, preset: str = "splice") -> list[str]:
    alignments = [a.split("\t") for a in align.to_sam(reference, query, preset=preset)]
    alignments_midsv = midsv.transform(alignments, midsv=False, cssplit=True, qscore=False)[0]
    if has_splice(alignments_midsv):
        alignments = [a.split("\t") for a in align.to_sam(reference, query, preset="map-ont")]
        alignments_midsv = midsv.transform(alignments, midsv=False, cssplit=True, qscore=False)[0]
    return alignments_midsv["CSSPLIT"].split(",")


def calculate_index_mapping(cssplits: list[str]) -> dict[int, int]:
    index_mapping = dict()
    query_index = 0
    range_inversion = set()
    for reference_index, element in enumerate(cssplits):
        # Inversion: If the last character is lower, it means we're in an inverted range
        if reference_index in range_inversion:
            continue
        if element[-1].islower():
            start = reference_index
            while element[-1].islower() and reference_index < len(cssplits):
                end = reference_index
                range_inversion.add(reference_index)
                reference_index += 1
                query_index += 1
                element = cssplits[reference_index]
            query_indexes = {i: j for i, j in zip(range(start, end + 1), range(end, start - 1, -1))}
            index_mapping.update(query_indexes)
        # Deletion or Unknown: Skip incrementing query_index
        elif element.startswith("-") or element == "N":
            continue
        # Substition: Just incrementing query_index
        elif element.startswith("*"):
            query_index += 1
        # Insertion: incrementing query_index by the number of insertions
        elif element.startswith("+"):
            query_index += element.count("+")
            index_mapping[reference_index] = query_index
            query_index += 1
        else:
            index_mapping[reference_index] = query_index
            query_index += 1
    return index_mapping


def get_index_mapping(TEMPDIR: str | Path) -> dict[str, dict[int, int]]:
    """
    Returns:
        dict(set): corresponding index of each allele based on control allele
        key: index of control allele
        value: index of other allele
    """
    path_fasta = [f for f in Path(TEMPDIR, "fasta").iterdir() if f.suffix != ".fai"]
    reference = select_control(path_fasta)
    INDEX_MAPPING = dict()
    for query in select_others(path_fasta):
        cssplits = get_cssplit(reference, query, preset="splice")
        allele = query.stem
        INDEX_MAPPING[allele] = calculate_index_mapping(cssplits)
    return INDEX_MAPPING


def save_index_mapping(TEMPDIR: str | Path) -> None:
    INDEX_MAPPING = get_index_mapping(TEMPDIR)
    with open(Path(TEMPDIR, "mutation_loci", "index_mapping.pickle"), "wb") as f:
        pickle.dump(INDEX_MAPPING, f)
