from __future__ import annotations

import re
import midsv

from pathlib import Path
from itertools import groupby
from DAJIN2.core import preprocess

###########################################################
# reverse complement to cssplits
###########################################################


def _reverse_cssplits(cssplits: list) -> list:
    for i, cs in enumerate(cssplits):
        if cs.startswith("+"):
            cssplits[i] = "+" + "|".join(cs.split("|")[::-1])
    return cssplits[::-1]


def _realign_insertion(cssplits: list) -> list:
    for i, cs in enumerate(cssplits):
        if not cs.startswith("+"):
            continue
        if re.search(rf"[{cs[1]}]", "[ACGTacgt]"):
            continue
        if i + 1 == len(cssplits):
            continue
        cs_current = cs.split("|")
        cssplits[i] = cs_current[0].replace("+", "")
        cssplits[i + 1] = "|".join(c[0] + c[-1] for c in cs_current[1:]) + "|" + cssplits[i + 1]
    return cssplits


def _complement_cssplit(cssplits: list) -> list:
    comp = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N", "a": "t", "c": "g", "g": "c", "t": "a", "n": "n"}
    for i, cs in enumerate(cssplits):
        op = cs[0]
        if op == "*":
            cssplits[i] = op + comp[cs[1]] + comp[cs[2]]
        elif op == "+":
            cssplits[i] = "|".join(c[0] + comp[c[-1]] for c in cs.split("|"))
        else:  # Match or Deletion or N
            cssplits[i] = op + comp[cs[-1]]
        cssplits[i] = cssplits[i].replace("NN", "N")
        cssplits[i] = cssplits[i].replace("nn", "n")
    return cssplits


def revcomp_cssplits(cssplits: list[str]) -> list[str]:
    cssplits_reversed = _reverse_cssplits(cssplits)
    cssplits_realigned = _realign_insertion(cssplits_reversed)
    cssplits_revcomped = _complement_cssplit(cssplits_realigned)
    return cssplits_revcomped


###########################################################
# group by mutation
###########################################################


def annotate_inversion(cssplits: list[str]) -> list[str]:
    return ["@" + cs if cs[-1].islower() else cs for cs in cssplits]


def group_by_mutation(cssplits: list[str]) -> list[list[str]]:
    return [list(group) for _, group in groupby(cssplits, key=lambda x: x[0])]


def flatten(lst: list[list]) -> list:
    results_flattend = []
    for result in lst:
        if isinstance(result[0], list):
            for res in result:
                results_flattend.append(res)
        else:
            results_flattend.append(result)
    return results_flattend


###########################################################
# report mutations
###########################################################


def _handle_match(group, genome, start, end, header, chromosome):
    end += len(group)
    start = end
    return None, start, end


def _handle_substitution(group, genome, start, end, header, chromosome):
    result = []
    for g in group:
        ref = g[1]
        mut = g[2]
        result.append([header, genome, chromosome, start, end, f"substitution: {ref}>{mut}"])
        end += 1
        start = end
    if len(result) == 1:
        result = result[0]
    return result, start, end


def _handle_deletion(group, genome, start, end, header, chromosome):
    end += len(group) - 1
    size = len(group)
    seq = "".join([g[-1] for g in group])
    result = [header, genome, chromosome, start, end, f"{size}bp deletion: {seq}"]
    end += 1
    start = end
    return result, start, end


def _handle_insertion(group, genome, start, end, header, chromosome):
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
    if len(result) == 1:
        result = result[0]
    end += 1
    start = end
    return result, start, end


def _handle_inversion(group, genome, start, end, header, chromosome):
    end += len(group) - 1
    size = len(group)
    seq = "".join([g[-1].upper() for g in group])
    result = [header, genome, chromosome, start, end, f"{size}bp inversion: {seq}"]
    end += 1
    start = end
    return result, start, end


def _handle_unknown(group, genome, start, end, header, chromosome):
    end += len(group) - 1
    size = len(group)
    result = [header, genome, chromosome, start, end, f"{size}bp unknown bases"]
    end += 1
    start = end
    return result, start, end


def report_mutations(cssplits_grouped, GENOME_COODINATES, header):
    genome = GENOME_COODINATES["genome"]
    chromosome = GENOME_COODINATES["chrom"]
    start = end = GENOME_COODINATES["start"]
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
                results.append(result)
    return flatten(results)


###########################################################
# main
###########################################################


def to_csv(TEMPDIR: Path | str, SAMPLE_NAME: str, GENOME_COODINATES: dict) -> None:
    results = [["Allele ID", "Genome", "Chromosome", "Start", "End", "Mutation"]]
    ref = Path(TEMPDIR, SAMPLE_NAME, "fasta", "control.fasta")
    for query in Path(TEMPDIR, "report", "FASTA", SAMPLE_NAME).iterdir():
        sam = preprocess.mapping.to_sam(ref, query)
        sam = [s.split("\t") for s in sam]
        midsv_sample = midsv.transform(sam, midsv=False, cssplit=True, qscore=False)[0]
        header = midsv_sample["QNAME"]
        cssplits = midsv_sample["CSSPLIT"].split(",")
        if GENOME_COODINATES.get("strand") == "-":
            cssplits = revcomp_cssplits(cssplits)
        cssplits_inversion = annotate_inversion(cssplits)
        cssplits_grouped = group_by_mutation(cssplits_inversion)
        result = report_mutations(cssplits_grouped, GENOME_COODINATES, header)
        results.extend(result)
    results_csv = "\n".join([",".join(map(str, r)) for r in results]) + "\n"
    path_output = Path(TEMPDIR, "report", "MUTATION_INFO", f"{SAMPLE_NAME}.csv")
    path_output.write_text(results_csv)
