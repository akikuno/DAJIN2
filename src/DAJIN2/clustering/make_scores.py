import re


def score_repeats(i, cssplit, scores, length, operant, inversion=False):
    j = i
    while j < length and cssplit[j].startswith(f"{operant}"):
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


def make_scores(sample_cssplit, diff_loci):
    diff_loci_cssplit = []
    for cssplit in sample_cssplit:
        cssplit = cssplit.split(",")
        diff_cs = []
        for locus in diff_loci:
            diff_cs.append(cssplit[locus])
        diff_loci_cssplit.append(diff_cs)
    length = len(diff_loci)
    for cssplit in diff_loci_cssplit:
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
        yield scores

