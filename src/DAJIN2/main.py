from __future__ import annotations
from itertools import groupby
import shutil
from pathlib import Path

# Custom modules
from src.DAJIN2.utils import cache_control
from src.DAJIN2.utils import io
from src.DAJIN2.utils import argparser

# from src.DAJIN2.classification import classification
# from src.DAJIN2.clustering import clustering
# from src.DAJIN2.consensus import consensus


#! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# def main():

########################################################################
# Argument parse
########################################################################

sample, control, allele, output, genome, debug, threads = argparser.parse()

# #* Point mutation
# sample, control, allele, output, genome, debug, threads = (
#     "examples/pm-tyr/barcode31.fq.gz",
#     "examples/pm-tyr/barcode32.fq.gz",
#     "examples/pm-tyr/design_tyr.fa",
#     "DAJIN_results",
#     "mm10",
#     True,
#     14
#     )


# #* 2-cut deletion
# sample, control, allele, output, genome, debug, threads = (
#     "examples/del-stx2/barcode25.fq.gz",
#     "examples/del-stx2/barcode30.fq.gz",
#     "examples/del-stx2/design_stx2.fa",
#     "DAJIN_results",
#     "mm10",
#     True,
#     14,
# )

# #* flox insertion
# sample, control, allele, output, genome, debug, threads = (
#     "examples/flox-cables2/AyabeTask1/barcode31.fq.gz",
#     "examples/flox-cables2/AyabeTask1/barcode42.fq.gz",
#     "examples/flox-cables2/AyabeTask1/design_cables2.fa",
#     "DAJIN_results",
#     "mm10",
#     True,
#     14,
# )

########################################################################
# Whether existing cached control
########################################################################

# CACHEDIR = os.path.join(tempfile.gettempdir(), "DAJIN")
# os.makedirs(CACHEDIR, exist_ok=True)

# if cache_control.exists(control, CACHEDIR):
#     IS_CACHED = True
# else:
#     cache_control.save_header(control, CACHEDIR)
#     IS_CACHED = False

###############################################################################
# Check and format inputs (sample/control/allele)
###############################################################################

from pathlib import Path
from importlib import reload
from src.DAJIN2.preprocess import format_input

reload(format_input)

# ------------------------------------------------------------------------------
# Check input path
# ------------------------------------------------------------------------------

if not Path(control).exists():
    raise FileNotFoundError(f"{control} is not found")
elif not Path(sample).exists():
    raise FileNotFoundError(f"{sample} is not found")
elif not Path(allele).exists():
    raise FileNotFoundError(f"{allele} is not found")

# ------------------------------------------------------------------------------
# Check formats (extensions and contents)
# ------------------------------------------------------------------------------

for file_name in sample, control:
    format_input.check_fastq_extension(file_name)

for file_name in sample, control:
    format_input.check_fastq_content(file_name)

format_input.check_fasta_content(allele)

# ------------------------------------------------------------------------------
# Formats
# ------------------------------------------------------------------------------

sample_name = format_input.extract_basename(sample)
control_name = format_input.extract_basename(control)
dict_allele = format_input.dictionize_allele(allele)

########################################################################
# Make directories
########################################################################

from pathlib import Path
import shutil

if Path(output).exists():
    shutil.rmtree(output)
Path(output).mkdir(exist_ok=True)

if Path(".tmpDAJIN").exists():
    shutil.rmtree(".tmpDAJIN")

for subdir in ["fasta", "fastq", "sam", "midsv"]:
    Path(".tmpDAJIN", subdir).mkdir(parents=True, exist_ok=True)

################################################################################
# Export fasta files as single-FASTA format
################################################################################

# TODO: use yeild, not export
from pathlib import Path

for header, sequence in dict_allele.items():
    contents = "\n".join([">" + header, sequence]) + "\n"
    output_fasta = Path(".tmpDAJIN", "fasta", f"{header}.fasta")
    output_fasta.write_text(contents)


###############################################################################
# Mapping with minimap2/mappy
###############################################################################

from pathlib import Path
from src.DAJIN2.preprocess import mappy_align
from importlib import reload

reload(mappy_align)

for input_fasta in Path(".tmpDAJIN", "fasta").glob("*.fasta"):
    fasta_name = input_fasta.name.replace(".fasta", "")
    for fastq, fastq_name in zip([control, sample], [control_name, sample_name]):
        # Todo: 並行処理で高速化！！
        SAM = mappy_align.to_sam(str(input_fasta), fastq)
        output_sam = Path(".tmpDAJIN", "sam", f"{fastq_name}_{fasta_name}.sam")
        output_sam.write_text("\n".join(SAM))

########################################################################
# MIDSV conversion
########################################################################

from pathlib import Path
import midsv

for sampath in Path(".tmpDAJIN", "sam").iterdir():
    output = Path(".tmpDAJIN", "midsv", f"{sampath.stem}.jsonl")
    sam = midsv.read_sam(sampath)
    midsv_jsonl = midsv.transform(sam, midsv=False, cssplit=True, qscore=False)
    midsv.write_jsonl(midsv_jsonl, output)


########################################################################
# Phread scoreが0.1以下のリードについて再分配する
########################################################################

# Under construction...

########################################################################
# Classify alleles
########################################################################

from src.DAJIN2 import classification
from importlib import reload

reload(classification)

classif_sample = classification.classify_alleles(sample_name)
classif_control = classification.classify_alleles(control_name)

########################################################################
# Detect Structural variants
########################################################################

from src.DAJIN2 import classification
from importlib import reload

reload(classification)

for classifs in [classif_sample, classif_control]:
    for classif in classifs:
        if classification.detect_sv(classif["CSSPLIT"], threshold=50):
            classif["SV"] = True
        else:
            classif["SV"] = False

classif_sample.sort(key=lambda x: (x["ALLELE"], x["SV"]))

# #? read check
# from collections import defaultdict

# d = defaultdict(int)
# for m in classif_sample:
#     d["total"] += 1
#     if m["SV"]:
#         d["SV"] += 1
#         # print(m["QNAME"], m["CSSPLIT"])
#     else:
#         d[m["ALLELE"]] += 1

# d

########################################################################
# Clustering
########################################################################

# -----------------------------------------------------------------------
# Extract significantly different base loci between Sample and Control
# -----------------------------------------------------------------------

import midsv
from pathlib import Path
from collections import defaultdict
from importlib import reload
from itertools import groupby
from src.DAJIN2 import clustering

reload(clustering)

dict_cssplit_control = defaultdict(list[dict])
for ALLELE in dict_allele.keys():
    path_control = Path(".tmpDAJIN", "midsv", f"{control_name}_{ALLELE}.jsonl")
    cssplit_control = [cs["CSSPLIT"] for cs in midsv.read_jsonl(path_control)]
    dict_cssplit_control[ALLELE] = cssplit_control

diffloci_by_alleles = defaultdict(list[dict])
for (ALLELE, SV), group in groupby(classif_sample, key=lambda x: (x["ALLELE"], x["SV"])):
    cssplit_sample = [record["CSSPLIT"] for record in group]
    cssplit_control = dict_cssplit_control[ALLELE]
    sequence = dict_allele[ALLELE]
    diffloci = clustering.screen_different_loci(cssplit_sample, cssplit_control, sequence, alpha=0.01, threshold=0.05)
    diffloci_by_alleles[f'{{"ALLELE": "{ALLELE}", "SV": {SV}}}'] = diffloci

# #? read check
# ALLELE = "flox"
# SV = False
# sequence = dict_allele[ALLELE]
# cssplit_sample = []
# for c in classif_sample:
#     if c["ALLELE"] == ALLELE and c["SV"] == SV:
#         cssplit_sample.append(c["CSSPLIT"])

# for s in cssplit_sample:
#     s.split(",")[1739]

# for c in classif_sample:
#     if c["QNAME"] == "880d5333-7a52-4dad-95fb-b778f06c1e5f":
#         cssplit = c["CSSPLIT"].split(",")
#         [i for i, cs in enumerate(cssplit) if cs.startswith("+A|+A|+C|+A|+T")]

# for c in classif_sample:
#     if c["QNAME"] == "6f536ce8-32d6-41d5-8b09-66bef425b9d8":
#         print(c)
#         cssplit = c["CSSPLIT"].split(",")
#         [i for i, cs in enumerate(cssplit) if cs.startswith("+A|+A|+C|+A|+T")]

# # diffloci_by_alleles
# for key, value in diffloci_by_alleles.items():
#     print(key, len(value))

# for key, value in diffloci_by_alleles.items():
#     if key == '{"ALLELE": "flox", "SV": False}':
#         for locus, score in value.items():
#             if any(abs(s) > 0.1 for s in score):
#                 print(locus, score)


# -----------------------------------------------------------------------
# Annotate scores to sample's reads
# -----------------------------------------------------------------------

from src.DAJIN2 import clustering
from collections import defaultdict

reload(clustering)

scores_by_alleles = defaultdict(list[dict])
for (ALLELE, SV), group in groupby(classif_sample, key=lambda x: (x["ALLELE"], x["SV"])):
    key = f'{{"ALLELE": "{ALLELE}", "SV": {SV}}}'
    cssplit_sample = [g["CSSPLIT"] for g in group]
    cssplit_control = dict_cssplit_control[ALLELE]
    diffloci = diffloci_by_alleles[key]
    scores_by_alleles[key] = clustering.make_scores(cssplit_sample, cssplit_control, diffloci)

# #? read check
# ALLELE = "control"
# SV = False
# cssplit_sample = []
# diffloci = diffloci_by_alleles[f'{{"ALLELE": "{ALLELE}", "SV": {SV}}}']
# for c in classif_sample:
#     if c["ALLELE"] == ALLELE and c["SV"] == SV:
#         cssplit_sample.append(c["CSSPLIT"])
#     if c["QNAME"] == "880d5333-7a52-4dad-95fb-b778f06c1e5f": # left-loxp
#         idx = len(cssplit_sample)-1

# x = list(clustering.make_scores(cssplit_sample, diffloci))
# x[idx]
# cssplit_sample[idx]

# -----------------------------------------------------------------------
# OPTICS clustering
# -----------------------------------------------------------------------

# ? read check
ALLELE = "control"
SV = False
cssplit_sample = []
cssplit_control = dict_cssplit_control[ALLELE]
qnames = []
diffloci = diffloci_by_alleles[f'{{"ALLELE": "{ALLELE}", "SV": {SV}}}']
for c in classif_sample:
    if c["ALLELE"] == ALLELE and c["SV"] == SV:
        cssplit_sample.append(c["CSSPLIT"])
        qnames.append(c["QNAME"])

d = defaultdict(int)
for c in classif_sample:
    ALLELE, SV = c["ALLELE"], c["SV"]
    key = f'{{"ALLELE": "{ALLELE}", "SV": {SV}}}'
    d[key] += 1


x = list(clustering.make_scores(cssplit_sample, cssplit_control, diffloci))

import numpy as np
import scipy.stats as st
from sklearn.cluster import Birch
from scipy.sparse import csr_matrix

# 3D to 2D
xarr = np.array(x)
xflatten = np.reshape(xarr, (len(xarr), -1))
xflatten = csr_matrix(xflatten)

# Clustering


def return_labels(xflatten, n_clusters=1):
    brc = Birch(n_clusters=n_clusters)
    brc.fit(xflatten)
    return brc.predict(xflatten)


def make_table(labels, n_clusters=1):
    count = defaultdict(int)
    for label in labels:
        count[label] += 1
    table = list(count.values())
    table.sort(reverse=True)
    if n_clusters == 1:
        table += [0]
    elif n_clusters == 2:
        pass
    else:
        table = table[:-1]
    return table


def chistatistic(s_table, c_table, threshold=0.05) -> float:
    chisq_value, _, ddof, _ = st.chi2_contingency([s_table, c_table])
    delta = sum(s_table + c_table) * threshold ** 2
    pval = 1 - st.ncx2.cdf(chisq_value, ddof, delta)
    return pval


prev_labels = return_labels(xflatten, n_clusters=1)
prev_table = make_table(prev_labels, n_clusters=1)

n_clusters = 2
while True:
    current_labels = return_labels(xflatten, n_clusters=n_clusters)
    current_table = make_table(current_labels, n_clusters=n_clusters)
    pval = chistatistic(prev_table, current_table, threshold=0.05)
    if pval > 0.05:
        labels = prev_labels
        break
    prev_table = current_table
    prev_labels = current_labels
    n_clusters += 1

# ? check

from collections import Counter

Counter(labels)

qa = []
for q, a in zip(qnames, labels.tolist()):
    qa.append({"QNAME": q, "LABEL": a})

qa.sort(key=lambda x: x["LABEL"])
from itertools import groupby

for LABEL, group in groupby(qa, key=lambda x: x["LABEL"]):
    with open(f"tmp_{LABEL}", "a") as f:
        f.write("^@\n")
        for g in group:
            # print(g["QNAME"])
            qname = g["QNAME"]
            f.write(f"{qname}\n")


# -----------------------------------------------------------------------
# (OLD) Score significantly different base loci between Sample and Control
# -----------------------------------------------------------------------

# import re
# from collections import defaultdict

# dict_diff_idsvn = defaultdict(list)
# classif_sample_groupby = groupby(classif_sample, key=lambda x: (x["ALLELE"], x["SV"]))
# for (ALLELE, SV), classif_sample_group in classif_sample_groupby:
#     diffloci = diffloci_by_alleles[f'{{"ALLELE": "{ALLELE}", "SV": {SV}}}']
#     diff_idsvn = [[0, 0, 0, 0, 0] for _ in diffloci]
#     for record in classif_sample_group:
#         cssplit_sample = record["CSSPLIT"].split(",")
#         for i, difflocus in enumerate(diffloci):
#             if cssplit_sample[difflocus].startswith("="):
#                 continue
#             if re.search(r"[acgtn]", cssplit_sample[difflocus]):
#                 diff_idsvn[i][3] += 1
#             elif cssplit_sample[difflocus].startswith("+"):
#                 diff_idsvn[i][0] += cssplit_sample[difflocus].count("+")
#             elif cssplit_sample[difflocus].startswith("-"):
#                 diff_idsvn[i][1] += 1
#             elif cssplit_sample[difflocus].startswith("*"):
#                 diff_idsvn[i][2] += 1
#             elif cssplit_sample[difflocus].startswith("N"):
#                 diff_idsvn[i][4] += 1
#         dict_diff_idsvn[f'{{"ALLELE": "{ALLELE}", "SV": {SV}}}'] = diff_idsvn

# dict_diff_idsvn_control = defaultdict(list)
# for ALLELE in dict_cssplit_control.keys():
#     group = dict_cssplit_control[ALLELE]
#     for SV in [True, False]:
#         diffloci = diffloci_by_alleles[f'{{"ALLELE": "{ALLELE}", "SV": {SV}}}']
#         diff_idsvn = [[0, 0, 0, 0, 0] for _ in diffloci]
#         for cssplit in group:
#             for i, difflocus in enumerate(diffloci):
#                 if cssplit[difflocus].startswith("="):
#                     continue
#                 if re.search(r"[acgtn]", cssplit[difflocus]):
#                     diff_idsvn[i][3] += 1
#                 elif cssplit[difflocus].startswith("+"):
#                     diff_idsvn[i][0] += cssplit[difflocus].count("+")
#                 elif cssplit[difflocus].startswith("-"):
#                     diff_idsvn[i][1] += 1
#                 elif cssplit[difflocus].startswith("*"):
#                     diff_idsvn[i][2] += 1
#                 elif cssplit[difflocus].startswith("N"):
#                     diff_idsvn[i][4] += 1
#             dict_diff_idsvn_control[f'{{"ALLELE": "{ALLELE}", "SV": {SV}}}'] = diff_idsvn

# count = 0
# for cssplit in group:
#     if cssplit[16].startswith("-"):
#         count += 1

# count

# dict_diff_idsvn_control.sort(key=lambda x: (x["ALLELE"], x["SV"]))


########################################################################
# レポート：アレル割合
########################################################################

########################################################################
# レポート：コンセンサス配列
########################################################################

########################################################################
# レポート：BAM
########################################################################

# ゲノム情報がなかったらリファレンスのFASTAにマップする
# ゲノム情報が入力されていたらその情報をもとにSAMの染色体位置情報を上書きする


if not debug:
    shutil.rmtree(".tmpDAJIN")

# if __name__ == "__main__":
#     main()
