from __future__ import annotations


def annotate_score(cssplits: list[list[str]], mutation_score: list[dict[str, float]]):
    scores = []
    for cssplit in cssplits:
        score = [0]
        for i in range(1, len(cssplit) - 1):
            if not mutation_score[i]:
                score.append(0)
                continue
            kmer = ",".join([cssplit[i - 1], cssplit[i], cssplit[i + 1]])
            score.append(mutation_score[i].get(kmer, 0))
        scores.append(score + [0])
    return scores

