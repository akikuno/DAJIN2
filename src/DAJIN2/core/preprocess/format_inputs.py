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
        # This is already a posix path, so just pass
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
# Fetch genome coodinate and chrom size
########################################################################


def fetch_coordinate(genome: str, ucsc_url: str, seq: str) -> dict:
    def _fetch_seq(seq: str):
        url_blat = f"{ucsc_url}/cgi-bin/hgBlat?db={genome}&type=BLAT&userSeq={seq}"
        request = urlopen(url_blat).read().decode("utf8").split("\n")
        matches = [x for x in request if "100.0%" in x]
        if not matches:
            raise AttributeError(f"{seq} is not found in {genome}")
        return matches[0].split()[-5:-1]

    seq_start, seq_end = seq[:1000], seq[-1000:]
    coordinate_start, coordinate_end = _fetch_seq(seq_start), _fetch_seq(seq_end)

    chromosome, strand = coordinate_start[0], coordinate_start[1]
    if strand == "+":
        start, end = int(coordinate_start[2]), int(coordinate_end[3])
    else:
        start, end = int(coordinate_end[2]), int(coordinate_start[3])

    return {"chr": chromosome, "start": start, "end": end, "strand": strand}


def fetch_chrom_size(chrom: str, genome: str, goldenpath_url: str) -> int:
    url = f"{goldenpath_url}/{genome}/bigZips/{genome}.chrom.sizes"
    request = urlopen(url).read().decode("utf8").split("\n")
    for req in request:
        req = req.split("\t")
        if chrom == req[0]:
            chrom_size = int(req[1])
    return chrom_size


def cache_coodinates_and_chromsize(TEMPDIR, GENOME, GENOME_COODINATES, CHROME_SIZE) -> None:
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
