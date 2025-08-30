from __future__ import annotations

import os
import ssl
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

import mappy
import pysam

from DAJIN2.utils.bed_handler import BEDError

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
# Check genome and UCSC server
# BAMファイルのゲノム座標の入手に必要であり、以下に、各Serverの役割を述べる：
## Blat: BAMファイル出力の際のゲノム座標の入手
## Goldenpath: chrom.sizeに必要
########################################################################


def fetch_html_without_verification(url: str, timeout: int = 10, retries: int = 3) -> str:
    """
    Fetch HTML with optional retries and clearer error messages.
    """
    context = ssl._create_unverified_context()  # Create an SSL context that temporarily disables verification

    for attempt in range(1, retries + 1):
        try:
            with urlopen(url, context=context, timeout=timeout) as response:
                return response.read().decode("utf-8", "ignore")
        except (HTTPError, URLError, TimeoutError) as e:
            if attempt < retries:
                time.sleep(1)
                continue
            else:
                raise TimeoutError(
                    f"Failed to fetch {url} after {retries} attempts (timeout={timeout}s). Last error: {e}"
                )


def format_url(key: str, url: str) -> str:
    return url + "/hg38" if key == "goldenpath" else url


def get_first_available_url(key: str, urls: list[str]) -> str | None:
    search_keys = {"gggenome": "hg38", "goldenpath": "bigZips"}
    return next(
        (url for url in urls if search_keys[key] in fetch_html_without_verification(format_url(key, url))), None
    )


def get_available_servers() -> dict[str, str]:
    server_lists = {
        "gggenome": [
            "https://gggenome.dbcls.jp/",
        ],
        "goldenpath": [
            "https://hgdownload.cse.ucsc.edu/goldenPath",
            "https://hgdownload.soe.ucsc.edu/goldenPath",
        ],
    }

    available_servers = {key: get_first_available_url(key, urls) for key, urls in server_lists.items()}

    error_messages = {
        "gggenome": "GGGenome servers are currently down. To avoid accessing the site, please consider specifying the -b/--bed option. Ref: 'https://github.com/akikuno/DAJIN2#using-bed-files-for-genomic-coordinates'",
        "goldenpath": "All UCSC GoldenPath servers are currently down. Please wait for a while and try again.",
    }

    for key, message in error_messages.items():
        if available_servers[key] is None:
            raise URLError(message)

    return available_servers


########################################################################
# BED file validation
########################################################################


def validate_bed_file(bed_path: str) -> None:
    """
    Validate BED file.

    Args:
        bed_path: Path to BED file
        genome: Optional genome assembly name

    Returns:
        Dictionary in DAJIN2 genome_coordinates format:
        {
            "genome": genome,
            "chrom": chromosome,
            "start": start_position,
            "end": end_position,
            "strand": strand,
            "chrom_size": chrom_size
        }

    Raises:
        FileNotFoundError: If BED file doesn't exist
        ValueError: If BED file is invalid
    """

    try:
        validate_file_existence(bed_path)

        # Check file extension
        path_bed = Path(bed_path)
        valid_extensions = [".bed", ".bed.gz"]
        file_suffixes = "".join(path_bed.suffixes).lower()

        if not any(file_suffixes.endswith(ext) for ext in valid_extensions):
            raise ValueError(f"BED file must have .bed or .bed.gz extension, got: {bed_path}")

        return None

    except BEDError as e:
        raise ValueError(f"Invalid BED file format: {e}")
    except Exception as e:
        raise ValueError(f"Error processing BED file {bed_path}: {e}")
