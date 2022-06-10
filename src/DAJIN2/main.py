import shutil
import sys
import os
import tempfile

# Custom modules
from src.DAJIN2.utils import cache_control
from src.DAJIN2.utils import io
from src.DAJIN2.utils import argparser
from src.DAJIN2.preprocess import mapping
from src.DAJIN2.preprocess import midsqc

# from src.DAJIN2.classification import classification
# from src.DAJIN2.clustering import clustering
# from src.DAJIN2.consensus import consensus


#! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# def main():

########################################################################
# Argument parse
########################################################################

sample, control, allele, output, genome, debug, threads = argparser.parse()

## Point mutation
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
#     "examples/nanosim/del-stx2/deletion.fq.gz",
#     "examples/nanosim/del-stx2/control.fq.gz",
#     "examples/nanosim/del-stx2/design_stx2.fa",
#     "DAJIN_results",
#     "mm10",
#     True,
#     14,
# )

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

########################################################################
# Make directories
########################################################################
import os

TMPDIR = ".tmpDAJIN"
SUBDIRS = ["fasta", "fastq", "sam", "midsconv", "midsqc"]
for sub in SUBDIRS:
    dir = os.path.join(TMPDIR, sub)
    os.makedirs(dir, exist_ok=True)

TMPDIR_PATHS = {dirname: os.path.join(TMPDIR, dirname) for dirname in os.listdir(TMPDIR)}

os.makedirs(output, exist_ok=True)

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
# Format inputs (sample/control/allele)
###############################################################################

import importlib
from src.DAJIN2.preprocess import format_input

importlib.reload(format_input)
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

# ------------------------------------------------------------------------------
# Export fasta files as single-FASTA format
# ------------------------------------------------------------------------------
# TODO: use yeild, not export
import pathlib

for header, sequence in dict_allele.items():
    contents = "\n".join([">" + header, sequence]) + "\n"
    output_fasta = pathlib.Path(".tmpDAJIN", "fasta", f"{header}.fasta")
    output_fasta.write_text(contents)


###############################################################################
# Mapping with minimap2/mappy
###############################################################################

import pathlib
from src.DAJIN2.preprocess import mappy_align

for input_fasta in pathlib.Path(".tmpDAJIN", "fasta").glob("*.fasta"):
    fasta_name = input_fasta.name.replace(".fasta", "")
    for fastq, fastq_name in zip([control, sample], [control_name, sample_name]):
        # Todo: 並行処理で高速化！！
        SAM = mappy_align.to_sam(str(input_fasta), fastq)
        output_sam = pathlib.Path(".tmpDAJIN", "sam", f"{fastq_name}_{fasta_name}.sam")
        output_sam.write_text("\n".join(SAM))

########################################################################
# MIDS conversion
########################################################################
# For development

import os
from src.DAJIN2.preprocess import midsconv
import importlib

importlib.reload(midsconv)

for samfile in os.listdir(".tmpDAJIN/sam"):
    sampath = os.path.join(".tmpDAJIN", "sam", samfile)
    output = os.path.join(".tmpDAJIN/midsconv", samfile.replace(".sam", ".csv"))
    midscsv = []
    for mids in midsconv.sam_to_mids(sampath, threads):
        # Extract full length reads
        x = midsconv.extract_full_length_reads(mids)
        if x:
            midscsv.append(x)
    with open(output, "w") as f:
        f.write("\n".join(midscsv))

########################################################################
# Phread scoreが0.1以下のリードについて再分配する
########################################################################

# Under construction...

########################################################################
# MIDS アレル分類・異常検知スコア
########################################################################

import os
import glob
from collections import defaultdict

mids_classification_score = defaultdict(list)

# midsfile = "barcode31_albino.csv"
for midsfile in glob.glob(".tmpDAJIN/midsconv/*"):
    tmp_name = os.path.basename(midsfile).replace(".csv", "").split("_")
    sample_name = "_".join(tmp_name[:-1])
    allele_type = tmp_name[-1]
    with open(midsfile) as f:
        midscsv = f.read().splitlines()
    read_len = len(midscsv[0].split(",")) - 1
    for mids in midscsv:
        match_score = 0
        read_id, *X = mids.split(",")
        for x in X:
            if x == "M":
                match_score += 1
            elif x[0].isdigit():
                match_score -= int(x[:-1])
            else:
                match_score -= 1
        match_score /= read_len
        mids_classification_score[sample_name, read_id].append([allele_type, match_score])

mids_classification_score2 = []
for key, val in mids_classification_score.items():
    val = sorted(val, key=lambda x: -x[1])[0]
    mids_classification_score2.append([*key, *val])

mids_classification_score2[-5:]
midscsv[0]
read_len = len(midscsv[0].split(",")) - 1
allele_type = "control"

for mids in midscsv:
    match_score = 0
    read_id, *X = mids.split(",")
    for x in X:
        if x == "M":
            match_score += 1
        elif x[0].isdigit():
            match_score -= int(x[:-1])
        else:
            match_score -= 1
    mids_classification_score[read_id].append([allele_type, match_score])

# mids_classification_score[read_id]
# = midscsv[0].split(",")[1:]

########################################################################
# MIDS アレル分類
########################################################################

########################################################################
# MIDS 異常検知
########################################################################


########################################################################
# MIDS クラスタリングスコア
########################################################################

from collections import defaultdict

midsscore = defaultdict(int)

read_num = len(midscsv)
for mids in midscsv:
    mids_split = mids.split(",")[1:]
    for i, x in enumerate(mids_split):
        midsscore[(i, x)] += 1 / read_num

for (i, mids), v in midsscore.items():
    if i == 828:
        print(mids, v)


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


if debug is False:
    shutil.rmtree(".tmpDAJIN")

#! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if __name__ == "__main__":
    main()
