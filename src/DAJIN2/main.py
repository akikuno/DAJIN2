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
# Whether existing cached control
########################################################################

# CACHEDIR = os.path.join(tempfile.gettempdir(), "DAJIN")
# os.makedirs(CACHEDIR, exist_ok=True)

# if cache_control.exists(control, CACHEDIR):
#     IS_CACHED = True
# else:
#     cache_control.save_header(control, CACHEDIR)
#     IS_CACHED = False

########################################################################
# Make directories
########################################################################

from pathlib import Path

Path(output).mkdir(exist_ok=True)

for subdir in ["fasta", "fastq", "sam", "midsv"]:
    Path(".tmpDAJIN", subdir).mkdir(parents=True, exist_ok=True)

###############################################################################
# Format inputs (sample/control/allele)
###############################################################################

from importlib import reload
from src.DAJIN2.preprocess import format_input

reload(format_input)
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
    output = Path(".tmpDAJIN", "midsv", f"{sampath.stem}.csv")
    sam = midsv.read_sam(sampath)
    midsv_jsonl = midsv.transform(sam)
    midsv.write_jsonl(midsv_jsonl, output)


########################################################################
# Phread scoreが0.1以下のリードについて再分配する
########################################################################

# Under construction...

########################################################################
# MIDSV scoring
########################################################################

from src.DAJIN2.preprocess import midsvscore
from importlib import reload

reload(midsvscore)
midsv_score_sample = midsvscore.extract_possible_allele_and_score(sample_name)
midsv_score_control = midsvscore.extract_possible_allele_and_score(control_name)

########################################################################
# Detect Structural variants
########################################################################

from src.DAJIN2.classification import detect_sv
from importlib import reload

reload(detect_sv)

for dict_midsvs in [midsv_score_sample, midsv_score_control]:
    for i, dict_midsv in enumerate(dict_midsvs):
        if detect_sv.is_sv(dict_midsv["CSSPLIT"]):
            dict_midsvs[i]["SV"] = True
        else:
            dict_midsvs[i]["SV"] = False


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
