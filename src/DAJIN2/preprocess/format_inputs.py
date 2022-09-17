from __future__ import annotations

from pathlib import Path
from urllib.request import urlopen
import re
import mappy

########################################################################
# Extract basename
########################################################################


def extract_basename(fastq_path: str) -> str:
    name = Path(fastq_path).name
    name = re.sub(r"\..*$", "", name)
    name = name.lstrip()
    if not name:
        raise AttributeError(f"{fastq_path} is not valid file name")
    else:
        return re.sub(r'[\\|/|:|?|.|,|\'|"|<|>|\| |]', "-", name)


########################################################################
# Convert allele file to dictionary type fasta format
########################################################################

# TODO TEST
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


########################################################################
# Fetch genome coodinate and chrom size
########################################################################

# TODO TEST
def fetch_coodinate(genome: str, ucsc_url: str, seq: str) -> dict:
    ucsc_blat = f"{ucsc_url}/cgi-bin/hgBlat?db={genome}&type=BLAT&userSeq={seq}"
    request = urlopen(ucsc_blat).read().decode("utf8").split("\n")
    coodinate = [x for x in request if x.count("100.0%")]
    if not coodinate:
        raise AttributeError(f"{seq} is not found in {genome}")
    else:
        coodinate = coodinate[0].split()
        return {"chr": coodinate[-5], "start": int(coodinate[-3]), "end": int(coodinate[-2]), "strand": coodinate[-4]}


def fetch_chrom_size(chrom: str, genome: str, goldenpath_url: str) -> int:
    url = f"{goldenpath_url}/{genome}/bigZips/{genome}.chrom.sizes"
    request = urlopen(url).read().decode("utf8").split("\n")
    for req in request:
        req = req.split("\t")
        if chrom == req[0]:
            chrom_size = int(req[1])
    return chrom_size
