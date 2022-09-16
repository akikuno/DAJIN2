from __future__ import annotations
from pathlib import Path
import re
import mappy
from urllib.request import urlopen
from urllib.error import URLError


def exists(input_file: str):
    if not Path(input_file).exists():
        raise FileNotFoundError(f"{input_file} is not found")


########################################################################
# Check if the sample is in the proper format.
########################################################################


def fastq_extension(fastq_path: str):
    correct_extension = False
    if re.search(r".fastq$|.fastq.gz$|.fq$|.fq.gz$", fastq_path):
        correct_extension = True
    if not correct_extension:
        raise AttributeError(f"{fastq_path} requires extensions either 'fastq', 'fastq.gz', 'fq' or 'fq.gz'")


# File content
def fastq_content(fastq_path: str):
    name, seq, qual = [], [], []
    for n, s, q in mappy.fastx_read(fastq_path):
        name.append(n)
        seq.append(s)
        qual.append(q)
    if not (len(name) == len(seq) == len(qual) > 0):
        raise AttributeError(f"{fastq_path} is not a FASTQ format")


def fasta_content(fasta_path: str):
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
# Check genome and UCSC server
########################################################################


def available_url(urls: list[str]) -> tuple[str, bool]:
    flag_fail = False
    url = ""
    for url in urls:
        try:
            _ = urlopen(url)
        except URLError:
            flag_fail = True
        else:
            flag_fail = False
            break
    return url, flag_fail


def available_genome(genome: str, ucsc_url: str):
    url = f"{ucsc_url}/cgi-bin/das/{genome}/dna?segment=1:1,10"
    response = urlopen(url)
    content = response.read()
    response.close()
    if not content:
        raise AttributeError(f"{genome} is not listed in UCSC genome browser")

