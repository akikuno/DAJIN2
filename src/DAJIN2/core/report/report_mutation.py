from __future__ import annotations

from pathlib import Path
from itertools import groupby


###########################################################
# Report mutation size in CSV format
###########################################################


def _revcomp_cssplits(cssplits: list[str]) -> list[str]:
    comp = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    for i, cs in enumerate(cssplits):
        op = cs[0]
        if op == "*":
            cssplits[i] = op + comp[cs[1]] + comp[cs[2]]
        elif op == "+":
            cssplits[i] = "|".join(c[0] + comp[c[-1]] for c in cs.split("|")[:-1])
        else:  # Match or Deletion or N
            cssplits[i] = op + comp[cs[-1]]
    return cssplits[::-1]


def _group_by_mutation(cssplits: list[str]) -> list[tuple[str]]:
    return [tuple(group) for _, group in groupby(cssplits, key=lambda x: x[0])]


def _report_mutations(cssplits_grouped, GENOME_COODINATES, header):
    chr_genome = GENOME_COODINATES["chr"]
    start = end = GENOME_COODINATES["start"]
    results = []

    def handle_match(group, start, end, header):
        end += len(group)
        start = end + 1
        return None, start, end

    def handle_substitution(group, start, end, header):
        end += 1
        ref = group[0][1]
        mut = group[0][2]
        result = [header, chr_genome, start, end, f"substitution: {ref}>{mut} "]
        start = end + 1
        return result, start, end

    def handle_deletion(group, start, end, header):
        end += len(group) - 1
        size = len(group)
        if size < 50:
            seq = "".join([g[-1] for g in group])
            result = [header, chr_genome, start, end, f"{size}bp deletion: {seq}"]
        else:
            result = [header, chr_genome, start, end, f"{size}bp deletion"]
        start = end + 1
        return result, start, end

    def handle_insertion(group, start, end, header):
        end += 1
        group = group[0]
        size = group.count("|")
        if size < 50:
            seq = "".join([g[-1] for g in group.split("|")[:-1]])
            result = [header, chr_genome, start, end, f"{size}bp insertion: {seq}"]
        else:
            result = [header, chr_genome, start, end, f"{size}bp insertion"]
        start = end + 1
        return result, start, end

    handlers = {"=": handle_match, "*": handle_substitution, "-": handle_deletion, "+": handle_insertion}

    for group in cssplits_grouped:
        for prefix, handler in handlers.items():
            if group[0].startswith(prefix):
                result, start, end = handler(group, start, end, header)
                if result:
                    results.append(result)
                break

    return results


def _to_csv(header: str, cons_per: list[dict], GENOME_COODINATES: dict) -> str:
    cons_cssplits = [max(cons, key=cons.get) for cons in cons_per]
    if GENOME_COODINATES["strand"] == "-":
        cons_cssplits = _revcomp_cssplits(cons_cssplits)
    cssplits_grouped = _group_by_mutation(cons_cssplits)
    return _report_mutations(cssplits_grouped, GENOME_COODINATES, header)


def to_csv(TEMPDIR: Path | str, SAMPLE_NAME: str, GENOME_COODINATES: dict, cons_percentage: dict) -> None:
    results = [["Allele ID", "Chromosome", "Start", "End", "Mutation"]]
    for header, cons_per in cons_percentage.items():
        result = _to_csv(header, cons_per, GENOME_COODINATES)
        results.extend(result)
    results_csv = "\n".join([",".join(map(str, r)) for r in results]) + "\n"
    path_output = Path(TEMPDIR, "report", "ALLELE_INFO", f"{SAMPLE_NAME}.csv")
    path_output.write_text(results_csv)
