from __future__ import annotations

import hashlib
import os
import ssl
import xml.etree.ElementTree as ET
from pathlib import Path
from urllib.error import URLError
from urllib.request import urlopen

import mappy
import pysam

########################################################################
# To accommodate cases where a user might input negative values or
# excessively large values, update the number of threads
# from "1" to "max cpu threads - 1".
########################################################################


def update_threads(threads: int) -> int:
    available_threads = os.cpu_count() - 1
    threads_updated = max(1, min(threads, available_threads))
    return threads_updated


########################################################################
# Check if the sample is in the proper format.
########################################################################


def validate_file_existence(path_file: str) -> None:
    if not Path(path_file).exists():
        raise FileNotFoundError(f"{path_file} is not found")


def return_file_extension(path_file: Path) -> str:
    extensions = [".fastq", ".fastq.gz", ".fq", ".fq.gz", ".fasta", ".fasta.gz", ".fa", ".fa.gz", ".bam"]
    file_suffixes = "".join(path_file.suffixes)  # Combine all suffixes into a single string

    for ext in extensions:
        if file_suffixes == ext:
            return ext

    raise ValueError(
        f"{path_file} requires extensions either .fastq, .fastq.gz, .fq, .fq.gz, .fasta, .fasta.gz, .fa, .fa.gz, or .bam"
    )


# Varidate if the file is in the proper format viewing top 100 lines
def validate_fastq_content(path_fastq: str) -> None:
    try:
        headers, seqs, quals = zip(*[(n, s, q) for i, (n, s, q) in enumerate(mappy.fastx_read(path_fastq)) if i < 100])
        # Remove empty elements
        headers = [header for header in headers if header]
        seqs = [seq for seq in seqs if seq]
        quals = [qual for qual in quals if qual]

        if not (len(headers) == len(seqs) == len(quals) > 0):
            raise ValueError

    except ValueError:
        raise ValueError(f"{path_fastq} is not a proper FASTQ format")


def validate_fasta_content(path_fasta: str, allele_file: bool = False) -> None:
    try:
        headers, seqs = zip(*[(n, s) for n, s, _ in mappy.fastx_read(path_fasta)])
        # Remove empty elements
        headers = [header for header in headers if header]
        seqs = [seq for seq in seqs if seq]

        if len(headers) != len(seqs) or not headers:
            raise ValueError

    except ValueError:
        raise ValueError(f"{path_fasta} is not a proper FASTA format")

    if len(headers) != len(set(headers)):
        raise ValueError(f"{path_fasta} must include unique identifiers")

    if len(seqs) != len(set(seqs)):
        raise ValueError(f"{path_fasta} must include unique DNA sequences")

    if allele_file and "control" not in headers:
        raise ValueError(f"One of the headers in the {path_fasta} must be '>control'")


def validate_bam_content(path_bam: str) -> None:
    try:
        _ = pysam.AlignmentFile(path_bam, "rb")
    except ValueError:
        raise ValueError(f"{path_bam} is not a proper BAM format")


def validate_files(SAMPLE: str, CONTROL: str, ALLELE: str) -> None:
    for path_directory in [CONTROL, SAMPLE]:
        extentions = {return_file_extension(path_fastx) for path_fastx in Path(path_directory).iterdir()}
        if len(extentions) == 1:
            extention = next(iter(extentions))
        else:
            raise ValueError(
                f"{path_directory} contains multiple extensions. Please check if there are any incorrect files."
            )

        fasta_extensions = {".fasta", ".fa", ".fasta.gz", ".fa.gz"}
        fastq_extentions = {".fastq", ".fq", ".fastq.gz", ".fq.gz"}
        if extention in fasta_extensions:
            _ = [validate_fasta_content(str(path_fastx)) for path_fastx in Path(path_directory).iterdir()]
        elif extention in fastq_extentions:
            _ = [validate_fastq_content(str(path_fastx)) for path_fastx in Path(path_directory).iterdir()]
        else:
            _ = [validate_bam_content(str(path_bam)) for path_bam in Path(path_directory).iterdir()]

    validate_file_existence(ALLELE)
    validate_fasta_content(ALLELE, allele_file=True)


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


def fetch_html_without_verification(url: str) -> str:
    context = ssl._create_unverified_context()  # Create an SSL context that temporarily disables verification
    with urlopen(url, context=context, timeout=10) as response:
        return response.read().decode("utf-8")


def format_url(key: str, url: str) -> str:
    return url + "/hg38" if key == "goldenpath" else url


def get_first_available_url(key: str, urls: list[str]) -> str | None:
    search_keys = {"blat": "BLAT Search Genome", "das": "GRCh38/hg38", "goldenpath": "bigZips"}
    return next(
        (url for url in urls if search_keys[key] in fetch_html_without_verification(format_url(key, url))), None
    )


def fetch_xml_without_verification(url: str) -> bytes:
    """Fetch XML data from a given URL."""
    context = ssl._create_unverified_context()  # Create an SSL context that temporarily disables verification
    with urlopen(url, context=context, timeout=10) as response:
        return response.read()


def extract_genome_ids_from_xml(xml_data: bytes) -> set[str]:
    """Extract genome IDs from XML data."""
    root = ET.fromstring(xml_data)
    return {cc.attrib["id"] for child in root for cc in child if cc.tag == "SOURCE"}


def get_genome_ids_in_ucsc(url_das: str) -> set[str]:
    """Get available genome IDs in UCSC."""
    xml_data = fetch_xml_without_verification(url_das)
    return extract_genome_ids_from_xml(xml_data)


def is_genome_id_available_in_ucsc(genome: str, url_das: str) -> bool:
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
            "https://hgdownload.soe.ucsc.edu/goldenPath",
        ],
    }

    available_servers = {key: get_first_available_url(key, urls) for key, urls in server_lists.items()}

    error_messages = {
        "blat": "All UCSC blat servers are currently down. Please wait for a while and try again.",
        "das": "All UCSC DAS servers are currently down. Please wait for a while and try again.",
        "goldenpath": "All UCSC GoldenPath servers are currently down. Please wait for a while and try again.",
    }

    for key, message in error_messages.items():
        if available_servers[key] is None:
            raise URLError(message)

    if not is_genome_id_available_in_ucsc(genome, available_servers["das"]):
        raise ValueError(f"{genome} is not listed. Available genomes are in {available_servers['das']}")

    del available_servers["das"]

    return available_servers
