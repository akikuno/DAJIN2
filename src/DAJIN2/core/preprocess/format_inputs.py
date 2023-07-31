from __future__ import annotations

import re
from pathlib import Path
from urllib.request import urlopen

import mappy
import wslPath

########################################################################
# Make directories
########################################################################


def make_directories(TEMPDIR: Path, SUBDIRS: list[str], SUBDIRS_REPORT: list[str], SAMPLE_NAME: str, CONTROL_NAME: str):
    for subdir in SUBDIRS:
        Path(TEMPDIR, subdir).mkdir(parents=True, exist_ok=True)
    for reportdir in SUBDIRS_REPORT:
        if reportdir == "MUTATION_INFO":
            Path(TEMPDIR, "report", reportdir).mkdir(parents=True, exist_ok=True)
        else:
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


def fetch_coordinate(GENOME_COODINATES: dict, GENOME_URLS: dict, seq: str) -> dict:
    def _fetch_seq(genome: str, seq: str, blat: str):
        url_blat = f"{blat}?db={genome}&type=BLAT&userSeq={seq}"
        request = urlopen(url_blat).read().decode("utf8").split("\n")
        matches = [x for x in request if "100.0%" in x]
        if not matches:
            raise AttributeError(f"{seq} is not found in {genome}")
        return matches[0].split()[-5:-1]

    genome = GENOME_COODINATES["genome"]
    blat = GENOME_URLS["blat"]
    seq_start, seq_end = seq[:1000], seq[-1000:]
    coordinate_start = _fetch_seq(genome, seq_start, blat)
    coordinate_end = _fetch_seq(genome, seq_end, blat)

    chromosome, strand = coordinate_start[0], coordinate_start[1]
    if strand == "+":
        start, end = int(coordinate_start[2]), int(coordinate_end[3])
    else:
        start, end = int(coordinate_end[2]), int(coordinate_start[3])

    GENOME_COODINATES["chr"] = chromosome
    GENOME_COODINATES["start"] = start
    GENOME_COODINATES["end"] = end
    GENOME_COODINATES["strand"] = strand
    return GENOME_COODINATES


def fetch_chrom_size(GENOME_COODINATES: dict, GENOME_URLS: dict) -> int:
    chrom = GENOME_COODINATES["chr"]
    genome = GENOME_COODINATES["genome"]
    url_goldenpath = GENOME_URLS["goldenpath"]
    url = f"{url_goldenpath}/{genome}/bigZips/{genome}.chrom.sizes"
    request = urlopen(url).read().decode("utf8").split("\n")
    for req in request:
        req = req.split("\t")
        if chrom == req[0]:
            chrom_size = int(req[1])
    GENOME_COODINATES["chrom_size"] = chrom_size
    return GENOME_COODINATES
