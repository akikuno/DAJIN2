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
# Extract significantly different base loci between Sample and Control
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
# Score significantly different base loci between Sample and Control
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

# -----------------------------------------------------------------------
# Annotate scores to sample's reads
# -----------------------------------------------------------------------

from copy import deepcopy
from collections import defaultdict

cluster_sample = deepcopy(classif_sample)

for c in cluster_sample:
    del c["CSSPLIT"]

classif_sample.sort(key=lambda x: (x["ALLELE"], x["SV"], x["QNAME"]))
cluster_sample.sort(key=lambda x: (x["ALLELE"], x["SV"], x["QNAME"]))
classif_groupby = groupby(classif_sample, key=lambda x: (x["ALLELE"], x["SV"]))
cluster_groupby = groupby(cluster_sample, key=lambda x: (x["ALLELE"], x["SV"]))

for ((ALLELE, SV), classif), (_, cluster) in zip(classif_groupby, cluster_groupby):
    keyname = f'{{"ALLELE": "{ALLELE}", "SV": {SV}}}'
    diffloci = allele_diffloci[keyname]
    diffscores = dict_diff_idsvn[keyname]
    for clas, clus in zip(classif, cluster):
        cs = clas["CSSPLIT"].split(",")
        cluster_score = []
        for difflocus, idsvnscore in zip(diffloci, diffscores):
            if cs[difflocus].startswith("="):
                cluster_score.append(0)
            elif re.search(r"[acgtn]", cs[difflocus]):
                cluster_score.append(idsvnscore[3])
            elif cs[difflocus].startswith("+"):
                cluster_score.append(idsvnscore[0])
            elif cs[difflocus].startswith("-"):
                cluster_score.append(idsvnscore[1])
            elif cs[difflocus].startswith("*"):
                cluster_score.append(idsvnscore[2])
            elif cs[difflocus].startswith("N"):
                cluster_score.append(idsvnscore[4])
        clus["SCORE"] = cluster_score

# #? read check
# for c in cluster_sample:
#     if c["QNAME"] == "6f536ce8-32d6-41d5-8b09-66bef425b9d8": # left-loxp
#         print(c)


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
