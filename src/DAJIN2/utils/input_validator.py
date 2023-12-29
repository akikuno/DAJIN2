from __future__ import annotations

import os
import re
import hashlib

from pathlib import Path
from urllib.error import URLError
from urllib.request import urlopen
import xml.etree.ElementTree as ET

import mappy


def update_threads(threads: int) -> int:
    available_threads = os.cpu_count() - 1
    threads_updated = max(1, min(threads, available_threads))
    return threads_updated


########################################################################
# Check if the sample is in the proper format.
########################################################################


def validate_file_existence(input_file: str):
    if not Path(input_file).exists():
        raise FileNotFoundError(f"{input_file} is not found")


def validate_fastq_extension(fastq_path: str):
    if not re.search(r".fastq$|.fastq.gz$|.fq$|.fq.gz$", fastq_path):
        raise ValueError(f"{fastq_path} requires extensions either 'fastq', 'fastq.gz', 'fq' or 'fq.gz'")


# Varidate if the file is in the proper format.
# See top 100 lines
def validate_fastq_content(fastq_path: str):
    try:
        names, seqs, quals = zip(*[(n, s, q) for i, (n, s, q) in enumerate(mappy.fastx_read(fastq_path)) if i < 100])
        if not (len(names) == len(seqs) == len(quals) > 0):
            raise ValueError
    except ValueError:
        raise ValueError(f"{fastq_path} is not a FASTQ format")


def validate_fasta_content(fasta_path: str):
    try:
        names, seqs = zip(*[(n, s) for n, s, _ in mappy.fastx_read(fasta_path)])
        if len(names) != len(seqs) or not names:
            raise ValueError
    except ValueError:
        raise ValueError(f"{fasta_path} is not a proper FASTA format")
    if len(names) != len(set(names)):
        raise ValueError(f"{fasta_path} must include unique identifiers")
    if len(seqs) != len(set(seqs)):
        raise ValueError(f"{fasta_path} must include unique DNA sequences")
    if "control" not in names:
        raise ValueError(f"One of the headers in the {fasta_path} must be '>control'")


def validate_files(SAMPLE: str, CONTROL: str, ALLELE: str) -> None:
    for path_fastx in [CONTROL, SAMPLE, ALLELE]:
        validate_file_existence(path_fastx)
    for path_fastx in [CONTROL, SAMPLE]:
        validate_fastq_extension(path_fastx)
        validate_fastq_content(path_fastx)
    validate_fasta_content(ALLELE)


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


def is_webpage_available(url: str) -> bool:
    try:
        with urlopen(url) as response:
            html = response.read().decode("utf-8")
            if "TITLE" not in html:
                return True
            title = next((h for h in html.split("\n") if "<TITLE>" in h), "")
            return "Error" not in title
    except URLError:
        return False


def get_first_available_url(urls: list[str]) -> str | None:
    return next((url for url in urls if is_webpage_available(url)), None)


def fetch_xml_data(url: str) -> bytes:
    """Fetch XML data from a given URL."""
    with urlopen(url) as response:
        return response.read()


def extract_genome_ids_from_xml(xml_data: bytes) -> set:
    """Extract genome IDs from XML data."""
    root = ET.fromstring(xml_data)
    return {cc.attrib["id"] for child in root for cc in child if cc.tag == "SOURCE"}


def get_genome_ids_in_ucsc(url_das: str) -> set:
    """Get available genome IDs in UCSC."""
    xml_data = fetch_xml_data(url_das)
    return extract_genome_ids_from_xml(xml_data)


def is_genome_in_ucsc_ids(genome: str, url_das: str) -> bool:
    genome_ids = get_genome_ids_in_ucsc(url_das)
    return genome in genome_ids


def validate_genome_and_fetch_urls(genome: str) -> dict[str, str]:
    server_lists = {
        "blat": [
            "https://genome.ucsc.edu/cgi-bin/hgBlat",
            "https://genome-asia.ucsc.edu/cgi-bin/hgBlat",
            "https://genome-euro.ucsc.edu/cgi-bin/hgBlat",
        ],
        "das": [
            "https://genome.ucsc.edu/cgi-bin/das/dsn/",
            "https://genome-asia.ucsc.edu/cgi-bin/das/dsn/",
            "https://genome-euro.ucsc.edu/cgi-bin/das/dsn",
        ],
        "goldenpath": [
            "https://hgdownload.cse.ucsc.edu/goldenPath",
            "http://hgdownload-euro.soe.ucsc.edu/goldenPath",
        ],
    }

    available_servers = {key: get_first_available_url(urls) for key, urls in server_lists.items()}

    error_messages = {
        "blat": "All UCSC blat servers are currently down. Please wait for a while and try again.",
        "das": "All UCSC DAS servers are currently down. Please wait for a while and try again.",
        "goldenpath": "All UCSC GoldenPath servers are currently down. Please wait for a while and try again.",
    }

    for key, message in error_messages.items():
        if not available_servers[key]:
            raise URLError(message)

    if not is_genome_in_ucsc_ids(genome, available_servers["das"]):
        raise ValueError(f"{genome} is not listed. Available genomes are in {available_servers['das']}")

    del available_servers["das"]

    return available_servers
