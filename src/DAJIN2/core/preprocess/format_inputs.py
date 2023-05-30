from __future__ import annotations

import os
import re
from pathlib import Path
from urllib.request import urlopen

import mappy
import midsv
import wslPath

########################################################################
# Make directories
########################################################################


def make_directories(TEMPDIR: Path, SUBDIRS: list[str], SUBDIRS_REPORT: list[str], SAMPLE_NAME: str, CONTROL_NAME: str):
    for subdir in SUBDIRS:
        Path(TEMPDIR, subdir).mkdir(parents=True, exist_ok=True)
    for reportdir in SUBDIRS_REPORT:
        Path(TEMPDIR, "report", reportdir, SAMPLE_NAME).mkdir(parents=True, exist_ok=True)
    Path(TEMPDIR, "report", "BAM", CONTROL_NAME).mkdir(parents=True, exist_ok=True)
    Path(TEMPDIR, "cache", ".igvjs").mkdir(parents=True, exist_ok=True)


########################################################################
# Convert Path
########################################################################


def convert_to_posix_path(path: str) -> str:
    try:
        path = wslPath.toPosix(path)
    except ValueError:
        pass
    return str(path)


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


def dictionize_allele(path_fasta: str) -> dict:
    header, sequence = [], []
    for name, seq, _ in mappy.fastx_read(path_fasta):
        name = name.lstrip()
        if not name:
            raise AttributeError(f"{path_fasta} contains an empty header")
        name = re.sub(r'[\\|/|:|?|.|,|\'|"|<|>|\| |]', "-", name)
        header.append(name)
        sequence.append(seq.upper())
    return {h: s for h, s in zip(header, sequence)}


########################################################################
# Update threads
########################################################################


def update_threads(threads):
    threads_updated = min(int(threads), os.cpu_count() - 1)
    threads_updated = max(1, threads_updated)
    return threads_updated


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


def cache_coodinates_and_chromsize(TEMPDIR, GENOME, GENOME_COODINATES, CHROME_SIZE):
    """
    Save (1) genome_symbol.txt, (2) genome_coodinates.jsonl, (3) chrome_size.txt
    """
    # Save info to the cache directory
    Path(TEMPDIR, "cache", "genome_symbol.txt").write_text(GENOME + "\n")
    midsv.write_jsonl([GENOME_COODINATES], Path(TEMPDIR, "cache", "genome_coodinates.jsonl"))
    Path(TEMPDIR, "cache", "chrome_size.txt").write_text(str(CHROME_SIZE))
    # Save info to the .igvjs directory
    Path(TEMPDIR, "report", ".igvjs", "genome_symbol.txt").write_text(GENOME + "\n")
    midsv.write_jsonl([GENOME_COODINATES], Path(TEMPDIR, "report", ".igvjs", "genome_coodinates.jsonl"))
