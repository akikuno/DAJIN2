import re
import numpy as np


def extract_cssplit_at_diffloci(cssplit_sample, diffloci):
    cssplit_diffloci = []
    for cssplit in cssplit_sample:
        cssplit = cssplit.split(",")
        diff_cs = []
        for locus in diffloci:
            diff_cs.append(cssplit[locus])
        cssplit_diffloci.append(diff_cs)
    return cssplit_diffloci


# Summation scores at each loci
def sum_scores(cssplit_diffloci):
    length = len(cssplit_diffloci[0])
    scores = [[0 for _ in range(7)] for _ in range(length)]
    for i, cssplit in enumerate(zip(*cssplit_diffloci)):
        for cs in cssplit:
            if cs.startswith("-"):
                scores[i][4] += 1
            elif cs.startswith("*"):
                scores[i][5] += 1
            elif cs.startswith("N"):
                scores[i][6] += 1
            elif cs.startswith("+"):
                cs_insertion = cs.split("|")
                if cs_insertion[-1].startswith("="):
                    scores[i][0] += 1
                elif cs_insertion[-1].startswith("-"):
                    scores[i][1] += 1
                elif cs_insertion[-1].startswith("*"):
                    scores[i][2] += 1
                elif cs_insertion[-1].startswith("N"):
                    scores[i][3] += 1
    return scores


def score_repeats(i, cssplit, scores, length, operant, inversion=False):
    j = i
    while j < length and cssplit[j].startswith(operant):
        j += 1
    score_repeat = j - i
    if inversion:
        score_repeat *= -1
    for idx in range(i, j):
        if operant == "-":
            scores[idx][4] = score_repeat
        if operant == "*":
            scores[idx][5] = score_repeat
        if operant == "N":
            scores[idx][6] = score_repeat
    return j - 1, scores


def score_repeats_insertion(i, cs, scores, inversion=False):
    cs_ins = cs.split("|")
    score = len(cs_ins)
    if inversion:
        score *= -1
    if cs_ins[-1].startswith("="):
        scores[i][0] = score
    elif cs_ins[-1].startswith("-"):
        scores[i][1] = score
    elif cs_ins[-1].startswith("*"):
        scores[i][2] = score
    else:
        scores[i][3] = score
    return i, scores


def make_scores(cssplit_sample, diffloci):
    cssplit_diffloci = extract_cssplit_at_diffloci(cssplit_sample, diffloci)
    scores_summed = sum_scores(cssplit_diffloci)
    scores_summed = np.array(scores_summed)
    length = len(diffloci)
    scores_all = []
    for cssplit in cssplit_diffloci:
        # IM, ID, IS, IN, D, S, N
        scores = [[0 for _ in range(7)] for _ in range(length)]
        i = 0
        while i < length:
            cs = cssplit[i]
            inversion = False
            if re.search(r"[acgtn]", cs):
                inversion = True
            if cs.startswith("-"):
                i, scores = score_repeats(i, cssplit, scores, length, "-", inversion)
            elif cs.startswith("*"):
                i, scores = score_repeats(i, cssplit, scores, length, "*", inversion)
            elif cs.startswith("N"):
                i, scores = score_repeats(i, cssplit, scores, length, "N", inversion)
            elif cs.startswith("+"):
                i, scores = score_repeats_insertion(i, cs, scores, inversion)
            i += 1
        scores = np.array(scores)
        scores = (scores * scores_summed).tolist()
        scores_all.append(scores)
    return scores_all
