from __future__ import annotations

import os
import re
import mappy

########################################################################
# Check if the sample is in the proper format.
########################################################################

# File extention: eithr 'fastq', 'fastq.gz', 'fq' or 'fq.gz'
def check_fastq_extension(fastq_path: str):
    correct_extension = False
    if re.search(r".fastq$|.fastq.gz$|.fq$|.fq.gz$", fastq_path):
        correct_extension = True
    if not correct_extension:
        raise AttributeError(f"{fastq_path} requires extensions either 'fastq', 'fastq.gz', 'fq' or 'fq.gz'")


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
        raise AttributeError(f"{fastq_path} is not a FASTQ format")


def check_fasta_content(fasta_path: str):
    name, seq = [], []
    for n, s, _ in mappy.fastx_read(fasta_path):
        name.append(n)
        seq.append(s)
    if not len(name) == len(seq) > 0:
        raise AttributeError(f"{fasta_path} is not a FASTA format")
    if len(name) > len(set(name)):
        raise AttributeError(f"{fasta_path} must include unique DNA sequences")
    if name.count("control") == 0:
        raise AttributeError(f"{fasta_path} must include a 'control' sequence")


########################################################################
# Extract basename
########################################################################


def extract_basename(fastq_path: str) -> str:
    name = os.path.basename(fastq_path)
    name = re.sub(r".fastq$|.fastq.gz$|.fq$|.fq.gz$", "", name)
    name = name.lstrip()
    if not name:
        raise AttributeError(f"{fastq_path} is an invalid file name")
    else:
        return re.sub(r'[\\|/|:|?|.|,|\'|"|<|>|\| |]', "-", name)


########################################################################
# Convert allele file to dictionary type fasta format
########################################################################


def dictionize_allele(allele_path: str) -> dict:
    header, sequence = [], []
    for name, seq, _ in mappy.fastx_read(allele_path):
        name = name.lstrip()
        if not name:
            raise AttributeError(f"{allele_path} contains an invalid header name")
        else:
            name = re.sub(r'[\\|/|:|?|.|,|\'|"|<|>|\| |]', "-", name)
        header.append(name)
        sequence.append(seq.upper())
    return {h: s for h, s in zip(header, sequence)}
