from __future__ import annotations

from itertools import groupby
from pathlib import Path

import pandas as pd

from DAJIN2.utils.cssplits_handler import revcomp_cssplits

###########################################################
# group by mutation
###########################################################


def annotate_inversion(midsv_tag: list[str]) -> list[str]:
    return ["@" + tag if tag.islower() else tag for tag in midsv_tag]


def group_by_mutation(midsv_tag: list[str]) -> list[list[str]]:
    return [list(group) for _, group in groupby(midsv_tag, key=lambda x: x[0])]


###########################################################
# Report mutations
###########################################################


def _handle_match(group, genome, start, end, header, chromosome) -> tuple[list[list[str]], int, int]:
    end += len(group)
    start = end
    return [None], start, end


def _handle_substitution(group, genome, start, end, header, chromosome) -> tuple[list[list[str]], int, int]:
    result = []
    for g in group:
        ref = g[1]
        mut = g[2]
        result.append([header, genome, chromosome, start, end, f"substitution: {ref}>{mut}"])
        end += 1
        start = end
    return result, start, end


def _handle_deletion(group, genome, start, end, header, chromosome) -> tuple[list[list[str]], int, int]:
    end += len(group) - 1
    size = len(group)
    seq = "".join([g[-1] for g in group])
    result = [header, genome, chromosome, start, end, f"{size}bp deletion: {seq}"]
    end += 1
    start = end
    return [result], start, end


def _handle_insertion(group, genome, start, end, header, chromosome) -> tuple[list[list[str]], int, int]:
    group = group[0]
    size = group.count("|")
    seq_insertion = "".join([g[-1] for g in group.split("|")[:-1]])
    seq_last = group.split("|")[-1]
    result = []
    result.append([header, genome, chromosome, start, end, f"{size}bp insertion: {seq_insertion}"])
    if seq_last.startswith("="):
        pass
    elif seq_last.startswith("-"):
        result.append([header, genome, chromosome, start, end, f"1bp deletion: {seq_last[-1]}"])
    end += 1
    start = end
    return result, start, end


def _handle_inversion(group, genome, start, end, header, chromosome) -> tuple[list[list[str]], int, int]:
    count_deletion = sum(1 for g in group if "-" in g)
    size = len(group) - count_deletion
    end += size - 1
    seq = "".join([g[-1].upper() for g in group if "-" not in g])
    result = [header, genome, chromosome, start, end, f"{size}bp inversion: {seq}"]
    end += 1
    start = end
    return [result], start, end


def _handle_unknown(group, genome, start, end, header, chromosome) -> tuple[list[list[str]], int, int]:
    end += len(group) - 1
    size = len(group)
    result = [header, genome, chromosome, start, end, f"{size}bp unknown bases"]
    end += 1
    start = end
    return [result], start, end


def report_mutations(cssplits_grouped: list[list[str]], GENOME_COORDINATES, header) -> list[list[str]]:
    genome = GENOME_COORDINATES["genome"]
    chromosome = GENOME_COORDINATES["chrom"]
    start = end = GENOME_COORDINATES["start"]
    handlers = {
        "=": _handle_match,
        "*": _handle_substitution,
        "-": _handle_deletion,
        "+": _handle_insertion,
        "@": _handle_inversion,
        "N": _handle_unknown,
    }
    results = []
    for group in cssplits_grouped:
        for prefix, handler in handlers.items():
            if group[0].startswith(prefix):
                result, start, end = handler(group, genome, start, end, header, chromosome)
                if prefix == "=":
                    continue
                results.extend(result)
    return [list(map(str, r)) for r in results]


###########################################################
# main
###########################################################


def export_to_csv(
    TEMPDIR: Path, SAMPLE_NAME: str, GENOME_COORDINATES: dict, cons_midsv_tags: dict[str, list[str]]
) -> None:
    results = []
    for header, cons_midsv_tag in cons_midsv_tags.items():
        if GENOME_COORDINATES.get("strand") == "-":
            cons_midsv_tag = revcomp_cssplits(cons_midsv_tag)
        cons_midsv_tag_inversion = annotate_inversion(cons_midsv_tag)
        cons_midsv_tag_grouped = group_by_mutation(cons_midsv_tag_inversion)
        result = report_mutations(cons_midsv_tag_grouped, GENOME_COORDINATES, header)
        results.extend(result)

    col_names = ["Allele ID", "Genome", "Chromosome", "Start", "End", "Mutation"]
    df_results = pd.DataFrame(results, columns=col_names).sort_values(by=["Allele ID", "Start"])

    path_output = Path(TEMPDIR, "report", "MUTATION_INFO", f"{SAMPLE_NAME}.csv")
    df_results.to_csv(path_output, index=False, encoding="utf-8", lineterminator="\n")
