import os
import re
import mappy
from src.DAJIN2.utils.exceptions import InputFileError
from src.DAJIN2.utils import io

########################################################################
# Check if the sample is in the proper format.
########################################################################

# File extention: eithr 'fastq', 'fastq.gz', 'fq' or 'fq.gz'
def check_fastq_extension(fastq_path: str):
    correct_extension = False
    if re.search(r".fastq$|.fastq.gz$|.fq$|.fq.gz$", fastq_path):
        correct_extension = True
    if not correct_extension:
        raise InputFileError(f"{fastq_path} requires extensions either 'fastq', 'fastq.gz', 'fq' or 'fq.gz'")


# File content
def check_fastq_content(fastq_path: str):
    name, seq, qual = [], [], []
    for n, s, q in mappy.fastx_read(fastq_path):
        name.append(n)
        seq.append(s)
        qual.append(q)
    if len(name) == len(seq) == len(qual) > 0:
        pass
    else:
        raise InputFileError(f"{fastq_path} is not a FASTQ format")


def check_fasta_content(fasta_path: str):
    name, seq = [], []
    for n, s, _ in mappy.fastx_read(fasta_path):
        name.append(n)
        seq.append(s)
    if not len(name) == len(seq) > 0:
        raise InputFileError(f"{fasta_path} is not a FASTA format")
    if len(name) > len(set(name)):
        raise InputFileError(f"{fasta_path} must include unique DNA sequences")
    if name.count("control") == 0:
        raise InputFileError(f"{fasta_path} must include a 'control' sequence")


########################################################################
# Extract basename
########################################################################


def extract_basename(fastq_path: str) -> str:
    name = os.path.basename(fastq_path)
    return re.sub(r".fastq$|.fastq.gz$|.fq$|.fq.gz$", "", name)


########################################################################
# Convert allele file to dictionary type fasta format
########################################################################


def dictionize_allele(allele: str) -> dict:
    header, sequence = [], []
    for name, seq, _ in mappy.fastx_read(allele):
        header.append(name)
        sequence.append(seq.upper())
    return {h: s for h, s in zip(header, sequence)}


########################################################################
# "TooLong" flag for excessively long arrays out of Fastq
########################################################################


# def annotate_TooLong_to_fastq(fastqpath: str, dict_allele: dict) -> None:
#     """
#     Annotate "TooLong" tag after a header when sequence lenth is more than 1.2 times of max sequence length
#     """
#     fastq = io.fread(fastqpath)
#     max_length = max(len(seq) for seq in dict_allele.values())
#     # Combine four lines into one list
#     fastq = [fastq[i : i + 4] for i in range(0, len(fastq), 4)]
#     fastq_anno = []
#     for f in fastq:
#         header = f[0].split()[0]
#         if len(f[1]) > max_length * 1.2:
#             header = "-".join([header, "TooLong"])
#         fastq_anno.append([header, *f[1:]])
#     return fastq_anno
