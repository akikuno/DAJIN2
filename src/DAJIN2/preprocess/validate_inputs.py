from __future__ import annotations

import hashlib
import re
from pathlib import Path
from urllib.error import URLError
from urllib.request import urlopen

import mappy


def _exist_file(input_file: str):
    if not Path(input_file).exists():
        raise FileNotFoundError(f"{input_file} is not found")


########################################################################
# Check if the sample is in the proper format.
########################################################################


def _fastq_extension(fastq_path: str):
    if not re.search(r".fastq$|.fastq.gz$|.fq$|.fq.gz$", fastq_path):
        raise AttributeError(f"{fastq_path} requires extensions either 'fastq', 'fastq.gz', 'fq' or 'fq.gz'")


# Varidate if the file is in the proper format. See top 100 lines
def _fastq_content(fastq_path: str):
    name, seq, qual = [], [], []
    for i, (n, s, q) in enumerate(mappy.fastx_read(fastq_path)):
        name.append(n)
        seq.append(s)
        qual.append(q)
        if i == 100:
            break
    if not (len(name) == len(seq) == len(qual) > 0):
        raise AttributeError(f"{fastq_path} is not a FASTQ format")


def _fasta_content(fasta_path: str):
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
        raise AttributeError(f"One of the headers in the {fasta_path} must be '>control'")


def validate_files(SAMPLE: str, CONTROL: str, ALLELE: str) -> None:
    _exist_file(CONTROL)
    _exist_file(SAMPLE)
    _exist_file(ALLELE)
    _fastq_extension(CONTROL)
    _fastq_content(CONTROL)
    _fastq_extension(SAMPLE)
    _fastq_content(SAMPLE)
    _fasta_content(ALLELE)


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


def _check_url_availabilities(urls: list[str]) -> list[bool]:
    availabilities = []
    for url in urls:
        try:
            _ = urlopen(url, timeout=10)
        except URLError:
            availabilities.append(False)
        else:
            availabilities.append(True)
    return availabilities


def _extract_available_urls(urls: list[str], availabilities: list[bool]) -> list[str]:
    available_urls = []
    for url, flag in zip(urls, availabilities):
        if flag:
            available_urls.append(url)
    return available_urls


def _is_genome_listed_in_UCSC(genome: str, ucsc_url: str) -> None:
    url = f"{ucsc_url}/cgi-bin/das/{genome}/dna?segment=1:1,10"
    response = urlopen(url, timeout=10)
    content = response.read()
    response.close()
    if not content:
        raise AttributeError(
            f"{genome} is not listed in UCSC genome browser. Available genomes are listed in {ucsc_url}/cgi-bin/das/dsn"
        )


def validate_genome_and_fetch_urls(GENOME: str) -> dict[str, str]:
    # Check UCSC Server is available
    URLS_UCSC = [
        "https://genome.ucsc.edu/",
        "https://genome-asia.ucsc.edu/",
        "https://genome-euro.ucsc.edu/",
    ]
    url_availabilities = _check_url_availabilities(URLS_UCSC)
    if not any(url_availabilities):
        raise URLError("All servers of UCSC Genome Browsers are currently down. Please wait for a while and try again.")
    URL_UCSC_AVAILABLE = _extract_available_urls(URLS_UCSC, url_availabilities)[0]
    # Check UCSC Download Server is available
    URLS_GOLDENPATH = [
        "https://hgdownload.cse.ucsc.edu/goldenPath",
        "http://hgdownload-euro.soe.ucsc.edu/goldenPath",
    ]
    url_availabilities = _check_url_availabilities(URLS_UCSC)
    if not any(url_availabilities):
        raise URLError("All servers of UCSC GoldenPath are currently down. Please wait for a while and try again.")
    URL_GOLDENPATH_AVAILABLE = _extract_available_urls(URLS_GOLDENPATH, url_availabilities)[0]
    # Check input genome is listed in UCSC
    _is_genome_listed_in_UCSC(GENOME, URL_UCSC_AVAILABLE)
    return {"ucsc": URL_UCSC_AVAILABLE, "goldenpath": URL_GOLDENPATH_AVAILABLE}
