from __future__ import annotations
from collections import Counter
import midsv
import numpy as np
import re
from pathlib import Path
from itertools import chain
from copy import deepcopy
import statsmodels.api as sm
import bisect
from collections import defaultdict
from sklearn.cluster import MeanShift
from scipy.spatial.distance import cosine
from itertools import groupby

# *Stx2 deletion: 2012-2739 (727 bases)
# tests/data/knockout/test_barcode25.fq.gz
# tests/data/knockout/test_barcode30.fq.gz
# tests/data/knockout/design_stx2.fa


def transpose(cssplits):
    return [list(cs) for cs in zip(*cssplits)]


def call_count(transpose_cssplits: list[list[str]]) -> list[dict[str:int]]:
    cssplit_counts = []
    for cssplit in transpose_cssplits:
        count = Counter(cssplit)
        count = dict(count)
        cssplit_counts.append(count)
    return cssplit_counts


def count_deletons_in_repaet(count_control, repeat_span):
    list_deletions = []
    count_deletions = []
    repeat_idx = 0
    repeat_start = repeat_span[repeat_idx][0]
    repeat_end = repeat_span[repeat_idx][1]
    for i, counts in enumerate(count_control):
        if repeat_start <= i < repeat_end:
            count_del = sum(val for key, val in counts.items() if key.startswith("-")) + 1
            count_deletions.append(count_del)
        elif i == repeat_end:
            list_deletions.append(count_deletions)
            repeat_idx += 1
            if repeat_idx == len(repeat_span):
                break
            repeat_start = repeat_span[repeat_idx][0]
            repeat_end = repeat_span[repeat_idx][1]
            count_deletions = []
            if repeat_start == i:
                count_del = sum(val for key, val in counts.items() if key.startswith("-"))
                count_deletions.append(count_del)
    return list_deletions


def find_repetitive_dels(count_control, count_sample, sequence) -> set[int]:
    repeat_regrex = r"A{4,}|C{4,}|G{4,}|T{4,}|N{4,}"
    repeat_span = list(x.span() for x in re.finditer(repeat_regrex, sequence))
    control_repdels = count_deletons_in_repaet(count_control, repeat_span)
    sample_repdels = count_deletons_in_repaet(count_sample, repeat_span)
    cossims = [1 - cosine(c, s) for c, s in zip(control_repdels, sample_repdels)]
    high_cossims = [list(range(span[0], span[1])) for cossim, span in zip(cossims, repeat_span) if cossim > 0.95]
    return set(chain.from_iterable(high_cossims))


def sampling(cnt: Counter, size: int):
    elements = []
    probs = []
    coverage = sum(cnt.values())
    for key, val in cnt.items():
        elements.append(key)
        probs.append(val / coverage)
    np.random.seed(1)
    samples = np.random.choice(a=elements, size=size, p=probs)
    return samples


def replace_repdels(transposed_cssplits: list[list[str]], repeat_dels: set):
    replased_cssplits = deepcopy(transposed_cssplits)
    for i, cssplits in enumerate(replased_cssplits):
        if i not in repeat_dels:
            continue
        cnt = Counter(cssplits)
        size = sum(1 for cs in cssplits if cs.startswith("-"))
        if size == 0:
            continue
        for key in list(cnt.keys()):
            if key.startswith("-"):
                del cnt[key]
        samples = sampling(cnt, size)
        iter_samples = iter(samples)
        for j, cs in enumerate(cssplits):
            if cs.startswith("-"):
                replased_cssplits[i][j] = next(iter_samples)
    return replased_cssplits


def replace_both_ends_n(transposed_cssplits: list[list[str]]):
    d_samples = defaultdict(iter)
    for i, cssplits in enumerate(transposed_cssplits):
        cnt = Counter(cssplits)
        if cnt["N"] == 0:  # No N
            continue
        if len(cnt) == 1 and cnt["N"] > 0:  # All N
            continue
        size = cnt["N"]
        del cnt["N"]
        samples = sampling(cnt, size)
        d_samples[i] = iter(samples)
    cssplits_replaced = [list(cs) for cs in zip(*transposed_cssplits)]
    for i, cssplits in enumerate(cssplits_replaced):
        for j, cs in enumerate(cssplits):
            if cs != "N":
                break
            cssplits_replaced[i][j] = next(d_samples[j])
    for i, cssplits in enumerate(cssplits_replaced):
        cssplits = cssplits[::-1]
        for j, cs in enumerate(cssplits):
            if cs != "N":
                break
            cssplits_replaced[i][len(cssplits) - 1 - j] = next(d_samples[len(cssplits) - 1 - j])
    cssplits_replaced = cssplits_replaced[::-1]
    cssplits_replaced = [list(cs) for cs in zip(*cssplits_replaced)]
    return cssplits_replaced


###############################################################################
# Discard common errors
###############################################################################


def call_percent(cssplit_counts: list[dict[str:int]]) -> list[dict[str:int]]:
    cssplit_percent = []
    coverage = sum(cssplit_counts[0].values())
    for counts in cssplit_counts:
        percent = {k: v / coverage * 100 for k, v in counts.items()}
        cssplit_percent.append(percent)
    return cssplit_percent


def subtract_percentage(percent_control, percent_sample) -> list[dict]:
    sample_subtracted = []
    for cont, samp in zip(percent_control, percent_sample):
        samp = Counter(samp)
        samp.subtract(Counter(cont))
        sample_subtracted.append(dict(samp))
    return sample_subtracted


def discard_common_error(sample_subtracted, threshold=0.5):
    sample_discarded = []
    for samp in sample_subtracted:
        remained = {k: v for k, v in samp.items() if v > threshold}
        sample_discarded.append(remained)
    return sample_discarded


def discard_match(sample_subtracted):
    sample_discarded = []
    for samp in sample_subtracted:
        remained = {k: v for k, v in samp.items() if not k.startswith("=")}
        sample_discarded.append(remained)
    return sample_discarded


###############################################################################
# annotate scores
###############################################################################


def annotate_scores(count_compensated, percent_discarded):
    scores = []
    for samp, mutation in zip(count_compensated, percent_discarded):
        score = []
        if not mutation:
            # # TODO 全インデックスにスコアを当てるとデバッグがし易い。本番では計算量低減のために以下の2行は削除する。
            # score = [0 for _ in range(len(samp))]
            # scores.append(score)
            continue
        for key, val in mutation.items():
            for s in samp:
                if s == key:
                    score.append(val)
                else:
                    score.append(0)
        scores.append(score)
    return transpose(scores)


###############################################################################
from src.DAJIN2.core import clustering
from importlib import reload

reload(clustering)

allele = "deletion"
sv = True
sequence = FASTA_ALLELES[allele]
# knockin_loci = KNOCKIN_LOCI[allele]

# Control
midsv_control = midsv.read_jsonl((Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_{allele}.jsonl")))
cssplits_control = [cs["CSSPLIT"].split(",") for cs in midsv_control]
# Sample
cssplits_sample = [cs["CSSPLIT"].split(",") for cs in classif_sample if cs["ALLELE"] == allele and cs["SV"] == sv]
cssplits_corrected_control, cssplits_corrected_sample = clustering.correct_cssplits.correct_cssplits(
    cssplits_control, cssplits_sample, sequence
)
mutation_score = clustering.make_score.make_score(cssplits_corrected_control, cssplits_corrected_sample)
scores_control = clustering.annotate_score.annotate_score(cssplits_corrected_control, mutation_score)
scores_sample = clustering.annotate_score.annotate_score(cssplits_corrected_sample, mutation_score)
scores = scores_sample + scores_control[:1000]
labels = clustering.return_labels.return_labels(scores, THREADS)
labels_control = labels[len(scores_sample) :]
labels_sample = labels[: len(scores_sample)]
labels_merged = clustering.merge_clusters.merge_clusters(labels_control, labels_sample)
Counter(labels_merged)


d = defaultdict(int)
for cs in cssplits_corrected_sample:
    x = sum(1 for c in cs if c == "N")
    d[x] += 1

sorted(d.items(), key=lambda x: x[0])
sorted(d.items(), key=lambda x: -x[1])

# labels = []
# max_label = 0
# for label in [[1, 1, 2, 2, 3], [1, 1, 1, 0, 0]]:
#     labels_reorder = clustering.reorder_labels(label, start=max_label)
#     max_label = max(labels_reorder)
#     labels += labels_reorder

# labels
# classif_sample.sort(key=lambda x: (x["ALLELE"], x["SV"]))
# tmp_classif_sample = []
# for (allele, sv), group in groupby(classif_sample, key=lambda x: (x["ALLELE"], x["SV"])):
#     count = 0
#     for g in list(group)[:10]:
#         tmp_classif_sample.append(g)

# cssplits_sample = [cs["CSSPLIT"].split(",") for cs in tmp_classif_sample if cs["ALLELE"] == allele and cs["SV"] == sv]

# x, xlabels = add_labels(classif_sample, CONTROL_NAME, FASTA_ALLELES, KNOCKIN_LOCI, TEMPDIR, THREADS)

# Counter(xlabels)
# d = defaultdict(int)
# Counter([xx["LABEL"] for xx in x])
# for xx in x:
#     d[xx["LABEL"]] += 1

# d

# original = [["A", "C", "N", "N"]]
# replaced = [["A", "C", "N", "N"]]

# for i, cssplits in enumerate(replaced):
#     cssplits_reverse = cssplits[::-1]
#     for j, cs in enumerate(cssplits_reverse):
#         if cs != "N":
#             break
#         print(len(cssplits_reverse) - 1 - j)

# Mask common mutations
# transpose_control = transpose(cssplits_corrected_control)
# transpose_sample = transpose(cssplits_corrected_sample)
# count_compensated_control = call_count(transpose_control)
# count_compensated_sample = call_count(transpose_sample)
# percent_control = call_percent(count_compensated_control)
# percent_sample = call_percent(count_compensated_sample)
# percent_subtraction = subtract_percentage(percent_control, percent_sample)
# percent_discarded = discard_common_error(percent_subtraction, 0.5)
# mutation_score = discard_match(percent_discarded)

# scores_control = annotate_scores(transpose_control, mutation_score)
# scores_sample = annotate_scores(transpose_sample, mutation_score)

# Clustering


# large del: 904d78ec7a0e
# medium del: 58a72e7c003a
# small del: 88b025c03911

"""
cat DAJINResults/.tempdir/test-stx2-deletion/midsv/test_barcode25_deletion.jsonl | grep 88b025c03911 | grep "-" # 259 dels
cat DAJINResults/.tempdir/test-stx2-deletion/midsv/test_barcode25_deletion.jsonl | grep 58a72e7c003a | grep "-" # 618 dels
minimap2 -ax map-ont tmp_cont.fa tmp_del.fa --cs
minimap2 -ax map-ont tmp_cont.fa tmp_inv.fa --cs=long  | wc -l
cat tests/data/knockout/design_stx2.fa | grep "deletion" -A1 > tmp_del.fa
cat tests/data/knockout/design_stx2.fa | grep "control" -A1 > tmp_cont.fa
minimap2 -t 30 -ax map-ont tmp_cont.fa tests/data/knockout/test_barcode25.fq.gz --cs=long > tmp_cont.sam
minimap2 -t 30 -ax map-ont tmp_del.fa tests/data/knockout/test_barcode25.fq.gz --cs=long > tmp_del.sam
minimap2 -t 30 -ax map-ont tmp_cont.fa tests/data/knockout/test_barcode25.fq.gz --splice  --cs=long > tmp_cont_splice.sam
minimap2 -t 30 -ax map-ont tmp_del.fa tests/data/knockout/test_barcode25.fq.gz  --splice --cs=long > tmp_del_splice.sam
cat tmp_cont.sam | grep -v "^@" | grep e0174c91bd7b | wc -l
cat tmp_del.sam | grep -v "^@" | grep e0174c91bd7b | wc -l
cat tmp_cont_splice.sam | grep -v "^@" | grep e0174c91bd7b | grep "~" | wc -l
cat tmp_del_splice.sam | grep -v "^@" | grep e0174c91bd7b | grep "~" | wc -l

cat tmp_del.sam | grep -e "^@" -e e0174c91bd7b | samtools sort > tmp_del_nosplice.bam && samtools index tmp_del_nosplice.bam
cat tmp_del_splice.sam | grep -e "^@" -e e0174c91bd7b | samtools sort > tmp_del_splice.bam && samtools index tmp_del_splice.bam

cat tmp_cont_splice.sam | grep -v "^@" | grep 88b025c03911 | awk '{print $(NF-1)}' | grep "~" | wc -l
cat tmp_del_splice.sam | grep -v "^@" | grep 88b025c03911 | awk '{print $(NF-1)}' | grep "~" | wc -l

cat tmp_cont.sam | grep -v "^@" | grep 4635f5e | wc -l
cat tmp_del.sam | grep -v "^@" | grep 4635f5e | wc -l
cat tmp_cont_splice.sam | grep -v "^@" | grep 4635f5e | grep "~"
cat tmp_del_splice.sam | grep -v "^@" | grep 4635f5e | grep "~"

cat tmp_cont_splice.sam | grep -v "^@" | grep ed618576505a | wc -l
cat tmp_del_splice.sam | grep -v "^@" | grep ed618576505a | wc -l

# delでリードが1つになってしまうものとそうでないものの違いはあるのか？
cat tmp_cont.sam | grep -v "^@" | cut -f 1 | sort | uniq -c | awk '$1 == 2 {print $2}' | sort > tmp_cont2
cat tmp_del.sam | grep -v "^@" | cut -f 1 | sort | uniq -c | awk '$1 == 2 {print $2}' | sort > tmp_del2
cat tmp_del.sam | grep -v "^@" | cut -f 1 | sort | uniq -c | awk '$1 == 1 {print $2}' | sort > tmp_del1

echo "^@" > tmp_cont2_del2
join tmp_cont2 tmp_del2 >> tmp_cont2_del2

echo "^@" > tmp_cont2_del1
join tmp_cont2 tmp_del1 >> tmp_cont2_del1

cat tmp_cont.sam | grep -f tmp_cont2_del1 | samtools sort > tmp_cont2_del1_nosplice_to_cont.bam && samtools index tmp_cont2_del1_nosplice_to_cont.bam
cat tmp_cont.sam | grep -f tmp_cont2_del2 | samtools sort > tmp_cont2_del2_nosplice_to_cont.bam && samtools index tmp_cont2_del2_nosplice_to_cont.bam
cat tmp_cont_splice.sam | grep -f tmp_cont2_del1 | samtools sort > tmp_cont2_del1_splice_to_cont.bam && samtools index tmp_cont2_del1_splice_to_cont.bam
cat tmp_cont_splice.sam | grep -f tmp_cont2_del2 | samtools sort > tmp_cont2_del2_splice_to_cont.bam && samtools index tmp_cont2_del2_splice_to_cont.bam


cat tmp_del.sam | grep -f tmp_cont2_del1 | samtools sort > tmp_cont2_del1_nosplice_to_del.bam && samtools index tmp_cont2_del1_nosplice_to_del.bam
cat tmp_del.sam | grep -f tmp_cont2_del2 | samtools sort > tmp_cont2_del2_nosplice_to_del.bam && samtools index tmp_cont2_del2_nosplice_to_del.bam
cat tmp_del_splice.sam | grep -f tmp_cont2_del1 | samtools sort > tmp_cont2_del1_splice_to_del.bam && samtools index tmp_cont2_del1_splice_to_del.bam
cat tmp_del_splice.sam | grep -f tmp_cont2_del2 | samtools sort > tmp_cont2_del2_splice_to_del.bam && samtools index tmp_cont2_del2_splice_to_del.bam
"""

# allele3 = "32c98e7b99b9"  # allele1に混ざったallele3
# allele4 = "65b9c02ba49a"  # 元来のallele1
# allele5 = "7fd963d0f30c"  # 元来のallele3

large_del = "904d78ec7a0e"
medium_del = "58a72e7c003a"
small_del = "88b025c03911"

allele = "deletion"
sv = True

xxx = [cs for cs in classif_sample if cs["ALLELE"] == allele and cs["SV"] == sv]

for i, classif in enumerate(xxx):
    if large_del in classif["QNAME"]:
        idx_large_del = i  # 863

for i, classif in enumerate(xxx):
    if medium_del in classif["QNAME"]:
        idx_medium_del = i  # 938

for i, classif in enumerate(xxx):
    if small_del in classif["QNAME"]:
        idx_small_del = i  # 938


xxx[idx_large_del]
xxx[idx_medium_del]
xxx[idx_small_del]
xxx[idx_large_del]["CSSPLIT"].split(",")[1314:1320]
xxx[idx_medium_del]["CSSPLIT"].split(",")[1314:1320]
cssplits_sample[idx_large_del][1314:1320]
cssplits_sample[idx_medium_del][1314:1320]

sorted(percent_sample[1320].items(), key=lambda x: -x[1])
sorted(percent_control[1320].items(), key=lambda x: -x[1])
sorted(percent_discarded[1320].items(), key=lambda x: -x[1])

mutation_score[1320]
i = 1317
cssplit = cssplits[i]
scores[idx_large_del][1314:1320]
scores[idx_medium_del][1314:1320]
scores_sample[idx_large_del][1314:1320]
scores_sample[idx_medium_del][1314:1320]
[round(s) for s in scores_sample[idx_large_del]]
[round(s) for s in scores_sample[idx_medium_del]]
[round(s) for s in scores_sample[idx_small_del]]
sum(scores_sample[idx_large_del])
sum(scores_sample[idx_medium_del])
sum(scores_sample[idx_small_del])

sum(1 for cs in cssplits_sample[idx_medium_del] if cs.startswith("-"))
sum(1 for cs in cssplits_sample[idx_small_del] if cs.startswith("-"))

sum([1 for i, s in enumerate(scores_sample[idx_medium_del]) if s > 10])
sum([1 for i, s in enumerate(scores_sample[idx_small_del]) if s > 10])

from collections import Counter

Counter(labels)
labels[idx_large_del]
labels[idx_medium_del]
labels[idx_small_del]

Counter(labels_merged)
