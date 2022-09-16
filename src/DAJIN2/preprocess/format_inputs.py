from __future__ import annotations

from pathlib import Path
import re
import mappy

########################################################################
# Extract basename
########################################################################


def extract_basename(fastq_path: str) -> str:
    name = Path(fastq_path).name
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
