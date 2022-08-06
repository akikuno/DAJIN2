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
# Output significantly different base loci between Sample and Control
# -----------------------------------------------------------------------

import midsv
from pathlib import Path
from collections import defaultdict
from importlib import reload
from itertools import groupby
from src.DAJIN2.clustering import find_difference

reload(find_difference)


classif_sample.sort(key=lambda x: (x["ALLELE"], x["SV"]))
classif_sample_groupby = groupby(classif_sample, key=lambda x: (x["ALLELE"], x["SV"]))

allele_diffloci = defaultdict(list)
for (ALLELE, SV), group in classif_sample_groupby:
    dict_control_cssplit = defaultdict(list[dict])
    if not dict_control_cssplit[ALLELE]:
        control_path = Path(".tmpDAJIN", "midsv", f"{control_name}_{ALLELE}.jsonl")
        control_cssplit = [cs["CSSPLIT"] for cs in midsv.read_jsonl(control_path)]
        dict_control_cssplit[ALLELE] = control_cssplit
    control_cssplit = dict_control_cssplit[ALLELE]
    sample_cssplit = []
    for record in group:
        sample_cssplit.append(record["CSSPLIT"])
    diffloci = find_difference.screen_different_loci(sample_cssplit, control_cssplit)
    allele_diffloci[f'{{"ALLELE": "{ALLELE}", "SV": {SV}}}'] = diffloci

# allele_diffloci


# -----------------------------------------------------------------------
# 変異部位のスコアリング
# -----------------------------------------------------------------------

import re
from collections import defaultdict

dict_diff_idsvn = defaultdict(list)
classif_sample_groupby = groupby(classif_sample, key=lambda x: (x["ALLELE"], x["SV"]))
for (ALLELE, SV), group in classif_sample_groupby:
    diffloci = allele_diffloci[f'{{"ALLELE": "{ALLELE}", "SV": {SV}}}']
    diff_idsvn = [[0, 0, 0, 0, 0] for _ in diffloci]
    for record in group:
        cs = record["CSSPLIT"].split(",")
        QNAME = record["QNAME"]
        for i, difflocus in enumerate(diffloci):
            if cs[difflocus].startswith("="):
                continue
            if re.search(r"[acgtn]", cs[difflocus]):
                diff_idsvn[i][3] += 1
            elif cs[difflocus].startswith("+"):
                diff_idsvn[i][0] += cs[difflocus].count("+")
            elif cs[difflocus].startswith("-"):
                diff_idsvn[i][1] += 1
            elif cs[difflocus].startswith("*"):
                diff_idsvn[i][2] += 1
            elif cs[difflocus].startswith("N"):
                diff_idsvn[i][4] += 1
        dict_diff_idsvn[f'{{"ALLELE": "{ALLELE}", "SV": {SV}}}'] = diff_idsvn


# locus = 2463
# locus = 1750
# ALLELE = "flox"
# SV = False
# total = 0
# counts = [0, 0]
# qnames = []
# for c in classif_sample:
#     if c["ALLELE"] == ALLELE and c["SV"] == SV:
#         total += 1
#         cssplit = c["CSSPLIT"].split(",")
#         if cssplit[locus].startswith("="):
#             counts[0] += 1
#         if cssplit[locus].startswith("-"):
#             counts[1] += 1
#             qnames.append(c["QNAME"])

# Path(f"tmp_qnames_{locus}_deletion.csv").write_text("^@\n" + "\n".join(qnames) + "\n")

# for c in classif_sample:
#     if c["QNAME"] == "6f536ce8-32d6-41d5-8b09-66bef425b9d8":
#         print(c)
# total
# counts

# dict_diff_idsvn[f'{{"ALLELE": "{ALLELE}", "SV": {SV}}}']

# for keys, values in dict_diff_idsvn.items():
#     ALLELE, SV, QNAME = eval(keys).values()
#     if QNAME == "af4aab64-d008-4b66-977e-d726d1245363":
#         print(values)

diffloci = allele_diffloci[f'{{"ALLELE": "{ALLELE}", "SV": {SV}}}']
diff_idsvn = [[0, 0, 0, 0, 0] for _ in range(len(diffloci))]
classif_sample_groupby = groupby(classif_sample, key=lambda x: (x["ALLELE"], x["SV"]))

"""
- controlとtargetで分散の違いをFisherの正確検定して、有意差のある塩基位置のみに注目する？
- Sampleでmatchが99.5%の塩基位置は解析対象から除く
- その後Kmerなどで計算量を削減する

read1:=C,+T|+G|=T,=A,=G,-C,=T,*GA,=G
read2:=C,+T|+G|=T,=A,=G,-C,=T,=A,=G
↓
# I, D, S, N
[
    {"read1": [[0,1,0,0,0,0,0,0],[0,0,0,0,1,0,0,0],[0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0]]},
    {"read2": [[0,1,0,0,0,0,0,0],[0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0]]}
]
↓
# persentage
[
    {"percentage": [[0,1,0,0,0,0,0,0],[0,0,0,0,1,0,0,0],[0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0]]},
    {"read2": [[0,1,0,0,0,0,0,0],[0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0]]}
]
"""

sample_cssplit = [m["CSSPLIT"] for m in midsv_score_sample]
control_cssplit = [m["CSSPLIT"] for m in midsv_score_control]

count_match = []
for i, (sample_cs, control_cs) in enumerate(zip(sample_cssplit, control_cssplit)):
    pass

from collections import defaultdict
from copy import deepcopy

allelereads = defaultdict(int)
midsscore = defaultdict(list)

for m in midsv_score_control:
    ALLELE = m["ALLELE"]
    CSSPLIT = m["CSSPLIT"].split(",")
    allelereads[ALLELE] += 1
    if not midsscore[ALLELE]:
        length = len(m["CSSPLIT"].split(","))
        midsscore[ALLELE] = [0] * length
    for i, cs in enumerate(CSSPLIT):
        if cs.startswith("="):
            midsscore[ALLELE][i] += 1

mids_match = deepcopy(midsscore)

for (allele, readnum), scores in zip(allelereads.items(), midsscore.values()):
    for i, score in enumerate(scores):
        mids_match[allele][i] = score / readnum * 100

print(*mids_match["control"][:100], sep="\n")

# read_num = len(midscsv)
# for mids in midscsv:
#     mids_split = mids.split(",")[1:]
#     for i, x in enumerate(mids_split):
#         midsscore[(i, x)] += 1 / read_num

# for (i, mids), v in midsscore.items():
#     if i == 828:
#         print(mids, v)


########################################################################
# クラスタリング
########################################################################

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
