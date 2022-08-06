from __future__ import annotations
import re
from copy import deepcopy
from itertools import groupby
from collections import defaultdict
import scipy.stats as st


def score_different_loci(sample_cssplit: list[str], control_cssplit: list[str]) -> list[dict]:
    """Scoring mutation (insertion, deletion, substitution, inversion, and unknow) at statistically significant loci between sample and control

    Args:
        sample_cssplit (list[str]): List of sample's CSSPLITs
        control_cssplit (list[str]): List of control's CSSPLITs

    Returns:
        list[dict]: Scores at significantly different base loci
    """
    # Create empty templates
    length = len(control_cssplit[0].split(","))
    # MIDSVN table
    sample_table = [[1, 1, 1, 1, 1, 1] for _ in range(length)]
    control_table = [[1, 1, 1, 1, 1, 1] for _ in range(length)]
    # Calculate match and mismatch
    for sample_cs, control_cs in zip(sample_cssplit, control_cssplit):
        sample_cs = sample_cs.split(",")
        control_cs = control_cs.split(",")
        for i, cs in enumerate(sample_cs):
            if re.search(r"[acgtn]", cs):
                sample_table[i][4] += 1
            elif cs.startswith("="):
                sample_table[i][0] += 1
            elif cs.startswith("+"):
                sample_table[i][1] += cs.count("+")
            elif cs.startswith("-"):
                sample_table[i][2] += 1
            elif cs.startswith("*"):
                sample_table[i][3] += 1
            elif cs.startswith("N"):
                sample_table[i][5] += 1
        for i, cs in enumerate(control_cs):
            if re.search(r"[acgtn]", cs):
                control_table[i][4] += 1
            elif cs.startswith("="):
                control_table[i][0] += 1
            elif cs.startswith("+"):
                control_table[i][1] += cs.count("+")
            elif cs.startswith("-"):
                control_table[i][2] += 1
            elif cs.startswith("*"):
                control_table[i][3] += 1
            elif cs.startswith("N"):
                control_table[i][5] += 1
    # for sample_cs in sample_cssplit:
    #     sample_cs = sample_cs.split(",")
    #     for i, cs in enumerate(sample_cs):
    #         if re.search(r"[acgtn]", cs):
    #             sample_table[i][4] += 1
    #         elif cs.startswith("="):
    #             sample_table[i][0] += 1
    #         elif cs.startswith("+"):
    #             sample_table[i][1] += cs.count("+")
    #         elif cs.startswith("-"):
    #             sample_table[i][2] += 1
    #         elif cs.startswith("*"):
    #             sample_table[i][3] += 1
    #         elif cs.startswith("N"):
    #             sample_table[i][5] += 1
    # for control_cs in control_cssplit:
    #     control_cs = control_cs.split(",")
    #     for i, cs in enumerate(control_cs):
    #         if re.search(r"[acgtn]", cs):
    #             control_table[i][4] += 1
    #         elif cs.startswith("="):
    #             control_table[i][0] += 1
    #         elif cs.startswith("+"):
    #             control_table[i][1] += cs.count("+")
    #         elif cs.startswith("-"):
    #             control_table[i][2] += 1
    #         elif cs.startswith("*"):
    #             control_table[i][3] += 1
    #         elif cs.startswith("N"):
    #             control_table[i][5] += 1
    # Calculate match and mismatch
    different_loci = []
    for i, (s, c) in enumerate(zip(sample_table, control_table)):
        # odds = (s[1] * c[0]) / (s[0] * c[1])
        pval = st.chi2_contingency([s, c])[1]
        if pval < 0.01:  # and odds > 1
            different_loci.append(i)
    n_sample = len(sample_cssplit)
    n_control = len(control_cssplit)
    diffloci_scores = defaultdict(list)
    for i in different_loci:
        sample_score = [score / n_sample for score in sample_table[i][1:]]
        control_score = [score / n_control for score in control_table[i][1:]]
        score = [s - c for s, c in zip(sample_score, control_score)]
        diffloci_scores[i] = score
    return diffloci_scores


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
            cluster_score = []
            for difflocus, score in diffloci_scores.items():
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
            clus["SCORE"] = cluster_score
    return cluster_sample

