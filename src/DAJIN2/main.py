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
# Check input path
###############################################################################

from pathlib import Path

if not Path(control).exists():
    raise FileNotFoundError(f"{control} is not found")
elif not Path(sample).exists():
    raise FileNotFoundError(f"{sample} is not found")
elif not Path(allele).exists():
    raise FileNotFoundError(f"{allele} is not found")


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
    output = Path(".tmpDAJIN", "midsv", f"{sampath.stem}.csv")
    sam = midsv.read_sam(sampath)
    midsv_jsonl = midsv.transform(sam)
    midsv.write_jsonl(midsv_jsonl, output)


########################################################################
# Phread scoreが0.1以下のリードについて再分配する
########################################################################

# Under construction...

########################################################################
# Classify alleles
########################################################################

from src.DAJIN2.classification import midsvscore
from importlib import reload

reload(midsvscore)

classif_sample = midsvscore.classify_alleles(sample_name)
classif_control = midsvscore.classify_alleles(control_name)

########################################################################
# Detect Structural variants
########################################################################

from src.DAJIN2.classification import detect_sv
from importlib import reload

reload(detect_sv)

for dict_midsvs in [classif_sample, classif_control]:
    for i, dict_midsv in enumerate(dict_midsvs):
        if detect_sv.is_sv(dict_midsv["CSSPLIT"]):
            dict_midsvs[i]["SV"] = True
        else:
            dict_midsvs[i]["SV"] = False

# #? read check

# p = Path("tmp_id.txt")
# p.write_text("^@\n")
# with p.open("a") as f:
#     for m in classif_sample:
#         f.write("^" + m["QNAME"] + "\n")

# for m in midsv_score_control:
#     if "7cdb4acdbb1b" in m["QNAME"]:
#         print(m)

# from collections import defaultdict

# d = defaultdict(int)
# for m in classif_control:
#     d["total"] += 1
#     if m["SV"]:
#         d["SV"] += 1
#         # print(m["QNAME"], m["CSSPLIT"])
#     else:
#         d[m["ALLELE"]] += 1

# d

########################################################################
# MIDS クラスタリングスコア
########################################################################

# -----------------------------------------------------------------------
# 変異部位のスコアリング
# -----------------------------------------------------------------------

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
