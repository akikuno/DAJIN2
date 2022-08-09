from __future__ import annotations
import re
from copy import deepcopy
from itertools import groupby
import scipy.stats as st


def make_table(sample_cssplit: list[str], control_cssplit: list[str]) -> tuple(list, list):
    # Create empty templates
    length = len(control_cssplit[0].split(","))
    # MIDSVN table
    sample_table = [[1, 1] for _ in range(length)]
    control_table = [[1, 1] for _ in range(length)]
    # Calculate match and mismatch
    for sample_cs, control_cs in zip(sample_cssplit, control_cssplit):
        sample_cs = sample_cs.split(",")
        control_cs = control_cs.split(",")
        for i, cs in enumerate(sample_cs):
            if cs.startswith("="):
                sample_table[i][0] += 1
            elif cs.startswith("+"):
                sample_table[i][1] += cs.count("+")
            else:
                sample_table[i][1] += 1
        for i, cs in enumerate(control_cs):
            if cs.startswith("="):
                control_table[i][0] += 1
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
    sample_cssplit: list[str], control_cssplit: list[str], sequence: str, alpha: float = 0.01, threshold: float = 0.05
) -> list[dict]:
    """Scoring mutation (insertion, deletion, substitution, inversion, and unknow) at statistically significant loci between sample and control

    Args:
        sample_cssplit (list[str]): List of sample's CSSPLITs
        control_cssplit (list[str]): List of control's CSSPLITs

    Returns:
        list[dict]: Scores at significantly different base loci
    """
    sample_table, control_table = make_table(sample_cssplit, control_cssplit)
    repeat_regrex = r"A{4,}|C{4,}|G{4,}|T{4,}|N{4,}"
    repeat_span = (x.span() for x in re.finditer(repeat_regrex, sequence))
    repeat_start, repeat_end = next(repeat_span)
    # Calculate match and mismatch
    different_loci = []
    s_repeat = [1, 1]
    c_repeat = [1, 1]
    for i, (s, c) in enumerate(zip(sample_table, control_table)):
        if i < repeat_start:
            pval = chistatistic(s, c, threshold)
            if pval < alpha:
                different_loci.append(i)
        elif repeat_start <= i <= repeat_end:
            s_repeat[0] += s[0]
            s_repeat[1] += s[1]
            c_repeat[0] += c[0]
            c_repeat[1] += c[1]
        elif i == repeat_end:
            pval = chistatistic(s_repeat, c_repeat, threshold)
            if pval < alpha:
                for j in range(repeat_start, repeat_end):
                    s, c = sample_table[j], control_table[j]
                    pval = chistatistic(s, c, threshold)
                    if pval < alpha:
                        different_loci.append(i)
            repeat_start, repeat_end = next(repeat_span)
            s_repeat = [1, 1]
            c_repeat = [1, 1]
        else:
            pval = chistatistic(s, c, threshold)
            if pval < alpha:
                different_loci.append(i)
    return different_loci


# s = sample_table[46]
# c = control_table[46]
# chisq_value, _, ddof, _ = st.chi2_contingency([s, c])
# delta = sum(s + c) * 0.05 ** 2
# pval = 1 - st.ncx2.cdf(chisq_value, ddof, delta)
# st.chisquare([s, c])

# s = sample_table[16:20]
# c = control_table[16:20]
# s_match = sum([ss[0] for ss in s])
# s_not = sum([sum(ss[1:]) for ss in s])

# c_match = sum([ss[0] for ss in c])
# c_not = sum([sum(ss[1:]) for ss in c])
# s = [s_match, s_not]
# c = [c_match, c_not]
# chisq_value, _, ddof, _ = st.chi2_contingency([s, c])
# delta = sum(s + c) * 0.1 ** 2
# pval = 1 - st.ncx2.cdf(chisq_value, ddof, delta)
# pval

# st.chi2_contingency([sample_table[16], control_table[16]])[1]
# st.chi2_contingency([sample_table[17], control_table[17]])
# control_table[16:20]
# sample_table[16:20]
# control_table[17]
# sample_table[17]
# # n_sample = len(sample_cssplit)
# # n_control = len(control_cssplit)
# diffloci_scores = defaultdict(list)
# for i in different_loci:
#     sample_score = [score / n_sample for score in sample_table[i][1:]]
#     control_score = [score / n_control for score in control_table[i][1:]]
#     score = [s - c for s, c in zip(sample_score, control_score)]
#     diffloci_scores[i] = score


def annotate_scores(classif_sample: list[dict], allele_diffloci: list[dict]) -> list[dict]:
    """Annotate scores to sample reads

    Args:
        classif_sample (list[dict]): 
        allele_diffloci (list[dict]): 

    Returns:
        list[dict]: Dist contains "SCORE"
    """
    classif_sample.sort(key=lambda x: (x["ALLELE"], x["SV"], x["QNAME"]))
    cluster_sample = deepcopy(classif_sample)
    for c in cluster_sample:
        del c["CSSPLIT"]
    classif_groupby = groupby(classif_sample, key=lambda x: (x["ALLELE"], x["SV"]))
    cluster_groupby = groupby(cluster_sample, key=lambda x: (x["ALLELE"], x["SV"]))
    for ((ALLELE, SV), classif), (_, cluster) in zip(classif_groupby, cluster_groupby):
        keyname = f'{{"ALLELE": "{ALLELE}", "SV": {SV}}}'
        diffloci_scores = allele_diffloci[keyname]
        for clas, clus in zip(classif, cluster):
            cs = clas["CSSPLIT"].split(",")
            cluster_cssplit = []
            cluster_score = []
            for difflocus, score in diffloci_scores.items():
                cluster_cssplit.append(cs[difflocus])
                if cs[difflocus].startswith("="):
                    cluster_score.append(0)
                elif re.search(r"[acgtn]", cs[difflocus]):
                    cluster_score.append(score[3])
                elif cs[difflocus].startswith("+"):
                    cluster_score.append(score[0])
                elif cs[difflocus].startswith("-"):
                    cluster_score.append(score[1])
                elif cs[difflocus].startswith("*"):
                    cluster_score.append(score[2])
                elif cs[difflocus].startswith("N"):
                    cluster_score.append(score[4])
            clus["CSSPLIT"] = cluster_cssplit
            clus["SCORE"] = cluster_score
    return cluster_sample

