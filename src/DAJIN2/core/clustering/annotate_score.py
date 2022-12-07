from __future__ import annotations


def annotate_score(cssplits: list[list[str]], mutation_score: list[dict[str, float]]):
    scores = []
    for cssplit in cssplits:
        score = []
        for i in range(1, len(cssplit) - 1):
            if not mutation_score[i]:
                continue
            kmer = ",".join([cssplit[i - 1], cssplit[i], cssplit[i + 1]])
            if kmer in mutation_score[i].keys():
                score.append(mutation_score[i][kmer])
            else:
                score.append(0)
        scores.append(score)
    return scores

