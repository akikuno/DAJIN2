from __future__ import annotations
import re
import scipy.stats as st


def replaceN(sample_cs: list[str]) -> list[str]:
    # Replace N to @ at the left ends
    for i, cs in enumerate(sample_cs):
        if cs != "N":
            break
        sample_cs[i] = "@"
    # Replace N to @ at the right ends
    sample_cs = sample_cs[::-1]
    for i, cs in enumerate(sample_cs):
        if cs != "N":
            break
        sample_cs[i] = "@"
    sample_cs = sample_cs[::-1]
    return sample_cs


def make_table(sample_cssplit: list[str], cssplit_control: list[str]) -> tuple(list, list):
    # Create empty templates
    length = len(cssplit_control[0].split(","))
    # MIDSVN table
    sample_table = [[1, 1] for _ in range(length)]
    control_table = [[1, 1] for _ in range(length)]
    # Calculate match and mismatch
    for sample_cs, control_cs in zip(sample_cssplit, cssplit_control):
        sample_cs = sample_cs.split(",")
        sample_cs = replaceN(sample_cs)
        control_cs = control_cs.split(",")
        control_cs = replaceN(control_cs)
        for i, cs in enumerate(sample_cs):
            if cs.startswith("="):
                sample_table[i][0] += 1
            elif cs == "@":
                sample_table[i][0] += 1
            elif cs.startswith("+"):
                sample_table[i][1] += cs.count("+")
            else:
                sample_table[i][1] += 1
        for i, cs in enumerate(control_cs):
            if cs.startswith("="):
                control_table[i][0] += 1
            elif cs == "@":
                sample_table[i][0] += 1
            elif cs.startswith("+"):
                control_table[i][1] += cs.count("+")
            else:
                control_table[i][1] += 1
    return sample_table, control_table


def chistatistic(s_table, c_table, threshold=0.05) -> float:
    chisq_value, _, ddof, _ = st.chi2_contingency([s_table, c_table])
    delta = sum(s_table + c_table) * threshold ** 2
    pval = 1 - st.ncx2.cdf(chisq_value, ddof, delta)
    return pval


def screen_different_loci(
    sample_cssplit: list[str], cssplit_control: list[str], sequence: str, alpha: float = 0.01, threshold: float = 0.05
) -> list[dict]:
    """Scoring mutation (insertion, deletion, substitution, inversion, and unknow) at statistically significant loci between sample and control

    Args:
        sample_cssplit (list[str]): List of sample's CSSPLITs
        cssplit_control (list[str]): List of control's CSSPLITs

    Returns:
        list[dict]: Scores at significantly different base loci
    """
    different_loci = []
    sample_table, control_table = make_table(sample_cssplit, cssplit_control)
    repeat_regrex = r"A{4,}|C{4,}|G{4,}|T{4,}|N{4,}"
    repeat_span = (x.span() for x in re.finditer(repeat_regrex, sequence))
    try:
        repeat_start, repeat_end = next(repeat_span)
        repeat_end -= 1
    except StopIteration:
        different_loci = []
        for i, (s, c) in enumerate(zip(sample_table, control_table)):
            pval = chistatistic(s, c, threshold)
            if pval < alpha:
                different_loci.append(i)
        return different_loci
    # Calculate match and mismatch
    s_repeat = [0, 0]
    c_repeat = [0, 0]
    for i, (s, c) in enumerate(zip(sample_table, control_table)):
        if i < repeat_start:
            pval = chistatistic(s, c, threshold)
            if pval < alpha:
                different_loci.append(i)
        elif repeat_start <= i < repeat_end:
            s_repeat[0] += s[0]
            s_repeat[1] += s[1]
            c_repeat[0] += c[0]
            c_repeat[1] += c[1]
        elif i == repeat_end:
            s_repeat[0] += s[0]
            s_repeat[1] += s[1]
            c_repeat[0] += c[0]
            c_repeat[1] += c[1]
            pval = chistatistic(s_repeat, c_repeat, threshold)
            if pval < alpha:
                for j in range(repeat_start, repeat_end + 1):
                    s, c = sample_table[j], control_table[j]
                    pval = chistatistic(s, c, threshold)
                    if pval < alpha:
                        different_loci.append(j)
            try:
                repeat_start, repeat_end = next(repeat_span)
                repeat_end -= 1
                s_repeat = [0, 0]
                c_repeat = [0, 0]
            except StopIteration:
                pass
        else:
            pval = chistatistic(s, c, threshold)
            if pval < alpha:
                different_loci.append(i)
    return different_loci


# def annotate_scores(classif_sample: list[dict], allele_diffloci: list[dict]) -> list[dict]:
#     """Annotate scores to sample reads

#     Args:
#         classif_sample (list[dict]):
#         allele_diffloci (list[dict]):

#     Returns:
#         list[dict]: Dist contains "SCORE"
#     """
#     classif_sample.sort(key=lambda x: (x["ALLELE"], x["SV"], x["QNAME"]))
#     cluster_sample = deepcopy(classif_sample)
#     for c in cluster_sample:
#         del c["CSSPLIT"]
#     classif_groupby = groupby(classif_sample, key=lambda x: (x["ALLELE"], x["SV"]))
#     cluster_groupby = groupby(cluster_sample, key=lambda x: (x["ALLELE"], x["SV"]))
#     for ((ALLELE, SV), classif), (_, cluster) in zip(classif_groupby, cluster_groupby):
#         keyname = f'{{"ALLELE": "{ALLELE}", "SV": {SV}}}'
#         diffloci_scores = allele_diffloci[keyname]
#         for clas, clus in zip(classif, cluster):
#             cs = clas["CSSPLIT"].split(",")
#             cluster_cssplit = []
#             cluster_score = []
#             for difflocus, score in diffloci_scores.items():
#                 cluster_cssplit.append(cs[difflocus])
#                 if cs[difflocus].startswith("="):
#                     cluster_score.append(0)
#                 elif re.search(r"[acgtn]", cs[difflocus]):
#                     cluster_score.append(score[3])
#                 elif cs[difflocus].startswith("+"):
#                     cluster_score.append(score[0])
#                 elif cs[difflocus].startswith("-"):
#                     cluster_score.append(score[1])
#                 elif cs[difflocus].startswith("*"):
#                     cluster_score.append(score[2])
#                 elif cs[difflocus].startswith("N"):
#                     cluster_score.append(score[4])
#             clus["CSSPLIT"] = cluster_cssplit
#             clus["SCORE"] = cluster_score
#     return cluster_sample

