from __future__ import annotations
from itertools import groupby
import shutil
from pathlib import Path
from sys import dllhandle
from this import d
from DAJIN2.consensus.module_consensus import join_list_of_dict

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

subdirectoris = ["fasta", "fastq", "sam", "midsv", "bam", "reports"]
for subdir in subdirectoris:
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
from src.DAJIN2 import clustering

reload(clustering)

dict_cssplit_control = defaultdict(list[dict])
for ALLELE in dict_allele.keys():
    path_control = Path(".tmpDAJIN", "midsv", f"{control_name}_{ALLELE}.jsonl")
    cssplit_control = [cs["CSSPLIT"] for cs in midsv.read_jsonl(path_control)]
    dict_cssplit_control[ALLELE] = cssplit_control

classif_sample.sort(key=lambda x: (x["ALLELE"], x["SV"]))
diffloci_by_alleles = defaultdict(list[dict])
for (ALLELE, SV), group in groupby(classif_sample, key=lambda x: (x["ALLELE"], x["SV"])):
    cssplit_sample = [record["CSSPLIT"] for record in group]
    cssplit_control = dict_cssplit_control[ALLELE]
    sequence = dict_allele[ALLELE]
    diffloci = clustering.screen_different_loci(cssplit_sample, cssplit_control, sequence, alpha=0.01, threshold=0.05)
    diffloci_by_alleles[f'{{"ALLELE": "{ALLELE}", "SV": {SV}}}'] = diffloci

# -----------------------------------------------------------------------
# Clustering
# -----------------------------------------------------------------------

from src.DAJIN2 import clustering
from copy import deepcopy

labels = []
label_start = 1
for (ALLELE, SV), group in groupby(classif_sample, key=lambda x: (x["ALLELE"], x["SV"])):
    key = f'{{"ALLELE": "{ALLELE}", "SV": {SV}}}'
    cssplit_sample = [g["CSSPLIT"] for g in group]
    diffloci = diffloci_by_alleles[key]
    scores = list(clustering.make_scores(cssplit_sample, diffloci))
    labels += [label + label_start for label in clustering.clustering(scores).tolist()]
    label_start = len(set(labels)) + 1

clust_sample = deepcopy(classif_sample)
for clust, label in zip(clust_sample, labels):
    clust["LABEL"] = label
    del clust["CSSPLIT"]

clust_sample.sort(key=lambda x: (x["ALLELE"], x["SV"], x["LABEL"]))

########################################################################
# Consensus call
########################################################################

from src.DAJIN2.consensus import module_consensus as consensus
from collections import defaultdict
from importlib import reload

reload(consensus)

path = Path(".tmpDAJIN", "midsv", f"{control_name}_control.jsonl")
cssplit_control = midsv.read_jsonl(path)

path = Path(".tmpDAJIN", "midsv", f"{sample_name}_control.jsonl")
cssplit_sample = midsv.read_jsonl(path)
cssplit_sample = consensus.join_listdicts(clust_sample, cssplit_sample, key="QNAME")
cssplit_sample.sort(key=lambda x: (x["ALLELE"], x["SV"], x["LABEL"]))

cons_percentage = defaultdict(list)
cons_sequence = defaultdict(list)
for (ALLELE, SV, LABEL), cssplits in groupby(cssplit_sample, key=lambda x: (x["ALLELE"], x["SV"], x["LABEL"])):
    cons_per = consensus.call_percentage(list(cssplits), cssplit_control)
    cons_seq = consensus.call_sequence(cons_per)
    key = f'{{"ALLELE": "{ALLELE}", "SV": {SV}, "LABEL": {LABEL}}}'
    cons_percentage[key] = cons_per
    cons_sequence[key] = cons_seq

########################################################################
# Report：アレル割合
# sample, allele name, #read, %read
########################################################################
from collections import defaultdict
from src.DAJIN2.report import report_af
import warnings

warnings.simplefilter("ignore")
reload(report_af)

# ----------------------------------------------------------
# All data
# ----------------------------------------------------------

df_clust_sample = report_af.all_allele(clust_sample, sample_name)
df_clust_sample.to_csv(".tmpDAJIN/reports/read_classification_all.csv", index=False)
df_clust_sample.to_excel(".tmpDAJIN/reports/read_classification_all.xlsx", index=False)

# ----------------------------------------------------------
# Summary data
# ----------------------------------------------------------

df_allele_frequency = report_af.summary_allele(clust_sample, sample_name, cons_sequence, dict_allele)
df_allele_frequency.to_csv(".tmpDAJIN/reports/read_classification_summary.csv", index=False)
df_allele_frequency.to_excel(".tmpDAJIN/reports/read_classification_summary.xlsx", index=False)

# ----------------------------------------------------------
# Visualization
# ----------------------------------------------------------
g = report_af.plot(df_allele_frequency)
g.save(filename=".tmpDAJIN/reports/tmp_output.png", dpi=350)
g.save(filename=".tmpDAJIN/reports/tmp_output.pdf")

########################################################################
# Report：各種ファイルフォーマットに出力
########################################################################

# FASTA
# for key, value in cons_percentage.items():
#     cons_fasta = consensus.call_fasta(key, value)
#     Path(f"tmp_{key}.fasta").write_text(cons_fasta)

# VCF

# HTML


########################################################################
# Report：BAM
########################################################################
import pysam

# ゲノム情報がなかったらリファレンスのFASTAにマップする
# ゲノム情報が入力されていたらその情報をもとにSAMの染色体位置情報を上書きする
if genome:
    pass

for path_sam in Path(".tmpDAJIN", "sam").glob("*_control.sam"):
    name = path_sam.stem
    pysam.sort("-o", f".tmpDAJIN/bam/{name}.bam", str(path_sam))
    pysam.index(f".tmpDAJIN/bam/{name}.bam")

########################################################################
# Report：IGV.js
########################################################################

# IGV.js（10本のリードのみ表示）のために各サンプル50本程度のリードのみを抽出する

for path_sam in Path(".tmpDAJIN", "sam").glob("*_control.sam"):
    sam = midsv.read_sam(path_sam)
    pysam.sort("-o", f".tmpDAJIN/bam/{name}.bam", str(path_sam))
    pysam.index(f".tmpDAJIN/bam/{name}.bam")


########################################################################
# Finish call
########################################################################

if not debug:
    shutil.rmtree(".tmpDAJIN")

# if __name__ == "__main__":
#     main()
