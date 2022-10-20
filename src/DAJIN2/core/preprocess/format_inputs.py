from __future__ import annotations

from pathlib import Path
from urllib.request import urlopen
import re
import mappy
import midsv

########################################################################
# Make directories
########################################################################


def make_directories(TEMPDIR: Path):
    subdirectoris = ["cache", "fasta", "sam", "midsv", "report", "result"]
    for subdir in subdirectoris:
        Path(TEMPDIR, subdir).mkdir(parents=True, exist_ok=True)
    reportdirectories = ["HTML", "FASTA", "VCF", "BAM", ".igvjs"]
    for reportdir in reportdirectories:
        Path(TEMPDIR, "report", reportdir).mkdir(parents=True, exist_ok=True)


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


def get_coodinates_and_chromsize(
    TEMPDIR, GENOME, DICT_ALLELE, UCSC_URL, GOLDENPATH_URL, IS_CACHE_GENOME: bool
) -> tuple:
    path_genome_coodinates = Path(TEMPDIR, "cache", "genome_coodinates.jsonl")
    path_chrome_size = Path(TEMPDIR, "cache", "chrome_size.txt")
    if IS_CACHE_GENOME:
        GENOME_COODINATES = midsv.read_jsonl(path_genome_coodinates)[0]
        CHROME_SIZE = int(path_chrome_size.read_text())
    else:
        GENOME_COODINATES = fetch_coodinate(GENOME, UCSC_URL, DICT_ALLELE["control"])
        CHROME_SIZE = fetch_chrom_size(GENOME_COODINATES["chr"], GENOME, GOLDENPATH_URL)
    return GENOME_COODINATES, CHROME_SIZE


def cache_coodinates_and_chromsize(TEMPDIR, GENOME, GENOME_COODINATES, CHROME_SIZE):
    # Save info to the cache directory
    Path(TEMPDIR, "cache", "genome_symbol.txt").write_text(GENOME + "\n")
    midsv.write_jsonl([GENOME_COODINATES], Path(TEMPDIR, "cache", "genome_coodinates.jsonl"))
    Path(TEMPDIR, "cache", "chrome_size.txt").write_text(str(CHROME_SIZE))
    # Save info to the .igvjs directory
    Path(TEMPDIR, "report", ".igvjs", "genome_symbol.txt").write_text(GENOME + "\n")
    midsv.write_jsonl([GENOME_COODINATES], Path(TEMPDIR, "report", ".igvjs", "genome_coodinates.jsonl"))
