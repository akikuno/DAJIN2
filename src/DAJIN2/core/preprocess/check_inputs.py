from __future__ import annotations
from pathlib import Path
import re
import mappy
import hashlib
from urllib.request import urlopen
from urllib.error import URLError


def exists(input_file: str):
    if not Path(input_file).exists():
        raise FileNotFoundError(f"{input_file} is not found")


########################################################################
# Check if the sample is in the proper format.
########################################################################


def fastq_extension(fastq_path: str):
    if not re.search(r".fastq$|.fastq.gz$|.fq$|.fq.gz$", fastq_path):
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
        raise AttributeError(f"{fasta_path} must include unique identifiers")
    if len(seq) > len(set(seq)):
        raise AttributeError(f"{fasta_path} must include unique DNA sequences")
    if name.count("control") == 0:
        raise AttributeError(f"{fasta_path} must include a 'control' sequence")


def check_files(SAMPLE: str, CONTROL: str, ALLELE: str) -> None:
    exists(CONTROL)
    exists(SAMPLE)
    exists(ALLELE)
    fastq_extension(CONTROL)
    fastq_content(CONTROL)
    fastq_extension(SAMPLE)
    fastq_content(SAMPLE)
    fasta_content(ALLELE)


########################################################################
# Check Cache
########################################################################


def exists_cached_control(CONTROL: str, TEMPDIR: Path) -> bool:
    PATH_CACHE_HASH = Path(TEMPDIR, "cache", "control_hash.txt")
    EXISTS_CACHE_CONTROL = False
    if PATH_CACHE_HASH.exists():
        current_hash = hashlib.sha256(Path(CONTROL).read_bytes()).hexdigest()
        cashed_hash = PATH_CACHE_HASH.read_text()
        if current_hash == cashed_hash:
            EXISTS_CACHE_CONTROL = True
    return EXISTS_CACHE_CONTROL


def exists_cached_genome(GENOME: str, TEMPDIR: Path, EXISTS_CACHE_CONTROL: bool) -> bool:
    PATH_CACHE_GENOME = Path(TEMPDIR, "cache", "genome_symbol.txt")
    EXISTS_CACHE_GENOME = False
    if GENOME and EXISTS_CACHE_CONTROL and PATH_CACHE_GENOME.exists():
        cashed_genome = PATH_CACHE_GENOME.read_text()
        if GENOME == cashed_genome:
            EXISTS_CACHE_GENOME = True
    return EXISTS_CACHE_GENOME


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
        raise AttributeError(
            f"{genome} is not listed in UCSC genome browser. Available genomes are in {ucsc_url}/cgi-bin/das/dsn"
        )


def check_and_fetch_genome(GENOME: str) -> tuple(str, str):
    # Check UCSC Server
    UCSC_URLS = [
        "https://genome.ucsc.edu/",
        "https://genome-asia.ucsc.edu/",
        "https://genome-euro.ucsc.edu/",
    ]
    UCSC_URL, flag_fail = available_url(UCSC_URLS)
    if flag_fail:
        raise URLError("UCSC Servers are currently down")
    # Check UCSC Download Server
    GOLDENPATH_URLS = [
        "https://hgdownload.cse.ucsc.edu/goldenPath",
        "http://hgdownload-euro.soe.ucsc.edu/goldenPath",
    ]
    GOLDENPATH_URL, flag_fail = available_url(GOLDENPATH_URLS)
    if flag_fail:
        raise URLError("UCSC Download Servers are currently down")
    # Check input genome
    available_genome(GENOME, UCSC_URL)
    return UCSC_URL, GOLDENPATH_URL
