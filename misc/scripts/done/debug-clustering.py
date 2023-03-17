import numpy as np


def pearson_corr(x, y):
    x_diff = x - np.mean(x)
    y_diff = y - np.mean(y)
    return np.dot(x_diff, y_diff) / (np.sqrt(sum(x_diff ** 2)) * np.sqrt(sum(y_diff ** 2)))


chistatistic([10000, 10000], [10000, 20000])

idx = 24
table_sample[idx]
table_control[idx]
chistatistic([coverage, sum(table_sample[idx])], [coverage, sum(table_control[idx])])
chistatistic([coverage, table_sample[idx][0]], [coverage, table_control[idx][0]])
chistatistic([coverage, table_sample[idx][1]], [coverage, table_control[idx][1]])
repeat_start = 461
repeat_end = 466
table_sample[repeat_start:repeat_end]
table_control[repeat_start:repeat_end]
len(cssplit_sample)
len(cssplit_control)

x = [x[0] for x in table_sample[repeat_start:repeat_end]]
y = [x[0] for x in table_control[repeat_start:repeat_end]]

pearson_corr(x, y)

SAMPLE, CONTROL, ALLELE, NAME, GENOME, DEBUG, THREADS = (
    "misc/data/barcode35.fastq.gz",
    "misc/data/barcode21.fastq.gz",
    "examples/flox-cables2/AyabeTask1/design_cables2.fa",
    "debug-clustering",
    "mm10",
    True,
    14,
)
