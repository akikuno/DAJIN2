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


def exists_cached_control(control: str, tempdir: Path) -> bool:
    path_cache_hash = Path(tempdir, "cache", "control_hash.txt")
    if path_cache_hash.exists():
        current_hash = hashlib.sha256(Path(control).read_bytes()).hexdigest()
        cashed_hash = path_cache_hash.read_text()
        if current_hash == cashed_hash:
            return True
    return False


def exists_cached_genome(genome: str, tempdir: Path, exists_cache_control: bool) -> bool:
    path_cache_genome = Path(tempdir, "cache", "genome_symbol.txt")
    if genome and exists_cache_control and path_cache_genome.exists():
        cashed_genome = path_cache_genome.read_text()
        if genome == cashed_genome:
            return True
    return False


########################################################################
# Check genome and UCSC server
########################################################################


def _check_url_availability(url: str) -> bool:
    try:
        _ = urlopen(url, timeout=10)
        return True
    except (URLError, TimeoutError):
        return False


def get_first_available_url(urls: list[str]) -> str | None:
    return next((url for url in urls if _check_url_availability(url)), None)


def is_genome_listed_in_UCSC(genome: str, ucsc_url: str) -> bool:
    url = f"{ucsc_url}/cgi-bin/das/{genome}/dna?segment=1:1,10"
    try:
        response = urlopen(url, timeout=10)
        return bool(response.read())
    except (URLError, TimeoutError):
        return False


def validate_genome_and_fetch_urls(genome: str) -> dict[str, str]:
    ucsc_blat_servers = [
        "https://genome.ucsc.edu/cgi-bin/hgBlat",
        "https://genome-asia.ucsc.edu/cgi-bin/hgBlat",
        "https://genome-euro.ucsc.edu/cgi-bin/hgBlat",
    ]
    ucsc_das_servers = [
        "https://genome.ucsc.edu/cgi-bin/das/dsn/",
        "https://genome-asia.ucsc.edu/cgi-bin/das/dsn/",
        "https://genome-euro.ucsc.edu/cgi-bin/das/dsn",
    ]
    goldenpath_servers = [
        "https://hgdownload.cse.ucsc.edu/goldenPath",
        "http://hgdownload-euro.soe.ucsc.edu/goldenPath",
    ]
    available_servers = {
        "blat": get_first_available_url(ucsc_blat_servers),
        "das": get_first_available_url(ucsc_das_servers),
        "goldenpath": get_first_available_url(goldenpath_servers),
    }

    if not available_servers["blat"]:
        raise URLError("All UCSC blat servers are currently down. Please wait for a while and try again.")

    if not available_servers["goldenpath"]:
        raise URLError("All UCSC GoldenPath servers are currently down. Please wait for a while and try again.")

    if not available_servers["das"]:
        raise URLError("All UCSC DAS servers are currently down. Please wait for a while and try again.")

    if not is_genome_listed_in_UCSC(genome, available_servers["das"]):
        raise AttributeError(f"{genome} is not listed. Available genomes are in {available_servers['das']}")

    return available_servers
