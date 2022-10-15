from __future__ import annotations
import re
import scipy.stats as st


def replaceNtoMatch(sample_cs: list[str]) -> list[str]:
    # Replace N to @ at the left ends
    for i, cs in enumerate(sample_cs):
        if cs != "N":
            break
        sample_cs[i] = "="
    # Replace N to @ at the right ends
    sample_cs = sample_cs[::-1]
    for i, cs in enumerate(sample_cs):
        if cs != "N":
            break
        sample_cs[i] = "="
    sample_cs = sample_cs[::-1]
    return sample_cs


def make_table(cssplit_sample: list[str], cssplit_control: list[str]) -> tuple(list, list):
    # Create empty templates
    length = len(cssplit_control[0].split(","))
    # table of Deletion, Others
    table_sample = [[1, 1] for _ in range(length)]
    table_control = [[1, 1] for _ in range(length)]
    # Calculate match and mismatch
    for sample_cs, control_cs in zip(cssplit_sample, cssplit_control):
        sample_cs = sample_cs.split(",")
        sample_cs = replaceNtoMatch(sample_cs)
        control_cs = control_cs.split(",")
        control_cs = replaceNtoMatch(control_cs)
        for i, cs in enumerate(sample_cs):
            if cs.startswith("="):
                continue
            if cs.startswith("-"):
                table_sample[i][0] += 1
            elif cs.startswith("+"):
                table_sample[i][1] += cs.count("+")
            else:
                table_sample[i][1] += 1
        for i, cs in enumerate(control_cs):
            if cs.startswith("="):
                continue
            if cs.startswith("-"):
                table_control[i][0] += 1
            elif cs.startswith("+"):
                table_control[i][1] += cs.count("+")
            else:
                table_control[i][1] += 1
    return table_sample, table_control


def mask_table_control(table_control, masks_control):
    for i, is_mask in enumerate(masks_control):
        if is_mask:
            table_control[i] = [1, 1]
    return table_control


def chistatistic(s_table, c_table, threshold=0.05) -> float:
    chisq_value, _, ddof, _ = st.chi2_contingency([s_table, c_table])
    delta = sum(s_table + c_table) * threshold ** 2
    pval = 1 - st.ncx2.cdf(chisq_value, ddof, delta)
    return pval


def screen_different_loci(
    cssplit_sample: list[str],
    cssplit_control: list[str],
    sequence: str,
    masks_control: list[bool],
    alpha: float = 0.01,
    threshold: float = 0.05,
) -> list[dict]:
    """Scoring mutation (insertion, deletion, substitution, inversion, and unknow)
    at significant loci between sample and control

    Args:
        cssplit_sample (list[str]): List of sample's CSSPLITs
        cssplit_control (list[str]): List of control's CSSPLITs

    Returns:
        list[dict]: Scores at significantly different base loci
    """
    different_loci = []
    coverage = min(len(cssplit_sample), len(cssplit_control))
    table_sample, table_control = make_table(cssplit_sample, cssplit_control)
    table_control = mask_table_control(table_control, masks_control)
    repeat_regrex = r"A{4,}|C{4,}|G{4,}|T{4,}|N{4,}"
    repeat_span = (x.span() for x in re.finditer(repeat_regrex, sequence))
    try:
        repeat_start, repeat_end = next(repeat_span)
        repeat_end -= 1
    except StopIteration:
        # Compare each nucleotide locus when sequence does not have a repeat
        for i, (s, c) in enumerate(zip(table_sample, table_control)):
            pval = chistatistic([coverage, sum(s)], [coverage, sum(c)], threshold)
            if pval < alpha:
                different_loci.append(i)
        return different_loci
    # Sum of deletion and other mismatchs
    repeat_del = [1, 1]
    repeat_other = [1, 1]
    for i, (s, c) in enumerate(zip(table_sample, table_control)):
        if repeat_start <= i < repeat_end:
            repeat_del[0] += s[0]
            repeat_del[1] += c[0]
            repeat_other[0] += s[1]
            repeat_other[1] += c[1]
        elif i == repeat_end:
            repeat_del[0] += s[0]
            repeat_del[1] += c[0]
            repeat_other[0] += s[1]
            repeat_other[1] += c[1]
            # statistics to deletion
            pval = chistatistic([coverage, repeat_del[0]], [coverage, repeat_del[1]], threshold)
            if pval < pow(alpha, 3):
                for j in range(repeat_start, repeat_end + 1):
                    s, c = table_sample[j], table_control[j]
                    pval = chistatistic([coverage, s[0]], [coverage, c[0]], threshold)
                    if pval < pow(alpha, 3):
                        different_loci.append(j)
            # statistics to others
            pval = chistatistic([coverage, repeat_other[0]], [coverage, repeat_other[1]], threshold)
            if pval < alpha:
                for j in range(repeat_start, repeat_end + 1):
                    s, c = table_sample[j], table_control[j]
                    pval = chistatistic([coverage, s[1]], [coverage, c[1]], threshold)
                    if pval < alpha and j not in different_loci:
                        different_loci.append(j)
            repeat_del = [1, 1]
            repeat_other = [1, 1]
            try:
                repeat_start, repeat_end = next(repeat_span)
                repeat_end -= 1
            except StopIteration:
                pass
        else:
            pval = chistatistic([coverage, sum(s)], [coverage, sum(c)], threshold)
            if pval < alpha:
                different_loci.append(i)
    return different_loci
