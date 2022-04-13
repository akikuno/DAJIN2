import shutil
import sys
import os
import tempfile

# Custom modules
from src.DAJIN2.utils import cache_control
from src.DAJIN2.utils import io
from src.DAJIN2.utils import argparser
from src.DAJIN2.preprocess import format_input
from src.DAJIN2.preprocess import mapping
from src.DAJIN2.preprocess import midsqc

# from src.DAJIN2.classification import classification
# from src.DAJIN2.clustering import clustering
# from src.DAJIN2.consensus import consensus


# For development
# import importlib


#! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# def main():
########################################################################
# Setting
########################################################################

# Check dependencies
# for dependence in ["minimap2", "samtools"]:
#     if not shutil.which(dependence):
#         print(f"{dependence} is not found. Please install it.")
#         sys.exit(1)

TMPDIR = ".tmpDAJIN"
SUBDIRS = ["fasta", "fastq", "sam", "midsconv", "midsqc"]
for sub in SUBDIRS:
    dir = os.path.join(TMPDIR, sub)
    os.makedirs(dir, exist_ok=True)

TMPDIR_PATHS = {dirname: os.path.join(TMPDIR, dirname) for dirname in os.listdir(TMPDIR)}

########################################################################
# Argument parse
########################################################################

sample, control, allele, output, genome, debug, threads = argparser.parse()

# sample, control, allele, output, genome, debug, threads = (
#     "examples/pm-tyr/barcode31.fq.gz",
#     "examples/pm-tyr/barcode32.fq.gz",
#     "examples/pm-tyr/design_tyr.fa",
#     "DAJIN_results",
#     "mm10",
#     True,
#     14
#     )

os.makedirs(output, exist_ok=True)

########################################################################
# Whether existing cached control
########################################################################

CACHEDIR = os.path.join(tempfile.gettempdir(), "DAJIN")
os.makedirs(CACHEDIR, exist_ok=True)

if cache_control.exists(control, CACHEDIR):
    IS_CACHED = True
else:
    cache_control.save_header(control, CACHEDIR)
    IS_CACHED = False

########################################################################
# Format inputs (sample/control/allele)
########################################################################

dict_allele = format_input.dictionize_allele(allele)

# Export fasta files as single-FASTA format
for header, sequence in dict_allele.items():
    contents = "\n".join([">" + header, sequence]) + "\n"
    with open(f".tmpDAJIN/fasta/{header}.fasta", "w") as f:
        _ = f.write(contents)

# for fastqpath in [sample, control]:
#     fastq_anno = format_input.annotate_TooLong_to_fastq(fastqpath, dict_allele)
#     basename = os.path.basename(fastqpath)
#     io.fwrite(fastq_anno, f".tmpDAJIN/fastq/{basename}")

########################################################################
# Mapping with minimap2/mappy
########################################################################
# importlib.reload(mapping)

from src.DAJIN2 import mappy2sam
import cstag
import pathlib

sample = "examples/pm-tyr/barcode31.fq.gz"

p = pathlib.Path(".tmpDAJIN/fasta")
for input_fasta in p.glob("*.fasta"):
    input_fasta = str(input_fasta)
    fasta_name = input_fasta.split("/")[-1].split(".f")[0]
    sample_name = sample.split("/")[-1].split(".f")[0]
    output_sam = f".tmpDAJIN/sam/{sample_name}_{fasta_name}.sam"
    # Todo: 並行処理で高速化！！
    SAM = mappy2sam(input_fasta, sample)
    with open(output_sam, "w") as f:
        f.write("\n".join(SAM))


########################################################################
# Mask CS tag
########################################################################

# qname_pos_cslong = []
# for sam in SAM:
#     if sam[0] == "@":
#         continue
#     sam = sam.split("\t")
#     qname, pos, cigar, qseq, qual, cs = sam[0], sam[3], sam[5], sam[9], sam[10], sam[11]
#     cslong = cstag.lengthen(cs, cigar, qseq)
#     cslong_masked = cstag.mask(cslong, cigar, qual, 5)
#     qname_pos_cslong.append([qname, pos, cslong_masked, qual])


########################################################################
# MIDS conversion
########################################################################
# importlib.reload(midsconv)
import os
from src.DAJIN2.preprocess import midsconv

sampath = ".tmpDAJIN/sam/barcode31_albino.sam"
threads = 20
# SAMDIR = os.path.join("tests", "samTomids", "input")
# SAMFILES = [os.path.join(SAMDIR, _) for _ in os.listdir(SAMDIR)]

for samfile in os.listdir(".tmpDAJIN/sam"):
    sampath = os.path.join(".tmpDAJIN", "sam", samfile)
    output = os.path.join(".tmpDAJIN/midsconv", samfile.replace(".sam", ".csv"))
    midscsv = midsconv.sam_to_mids(sampath, threads)
    with open(output, "w") as f:
        f.write("\n".join(midscsv))

########################################################################
# MIDS QC filtering
# 完全長リードのみを取り出す：両端から50bp連続して"="であるリードを除く
# Phread scoreが0.1以下のリードについて再分配する
########################################################################

print(mids_csv[0])

if debug is False:
    shutil.rmtree(".tmpDAJIN")

#! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if __name__ == "__main__":
    main()


cstag = "cs:Z:=ACGT*ag=C-g=T+t=ACGT"

cstag = samdict["cstag"]

cstags = ["MMMM", "S", "M", "Dg", "M", "It", "MMMM"]
