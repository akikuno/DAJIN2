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


# from __future__ import annotations


# def transpose(cssplits: list[list[str]]) -> list[list[str]]:
#     return [list(cs) for cs in zip(*cssplits)]


# def annotate_score(cssplits: list[list[str]], mutation_score: list[dict[str, float]]):
#     transpose_cssplits = transpose(cssplits)
#     scores = []
#     for samp, mutation in zip(transpose_cssplits, mutation_score):
#         score = []
#         if not mutation:
#             continue
#         for key, val in mutation.items():
#             for s in samp:
#                 if s == key:
#                     score.append(val)
#                 else:
#                     score.append(0)
#         scores.append(score)
#     return transpose(scores)

