from __future__ import annotations

import ssl
import time
from urllib.error import HTTPError, URLError
from urllib.request import urlopen


def fetch_bed_without_verification(url: str, timeout: int = 10, retries: int = 3) -> list[str]:
    context = ssl._create_unverified_context()  # Create an SSL context that temporarily disables verification

    for attempt in range(1, retries + 1):
        try:
            with urlopen(url, context=context, timeout=timeout) as response:
                return [res.split("\t") for res in response.read().decode("utf-8").split("\n") if res]
        except (HTTPError, URLError, TimeoutError) as e:
            if attempt < retries:
                time.sleep(1)
                continue
            else:
                raise TimeoutError(
                    f"Failed to fetch {url} after {retries} attempts (timeout={timeout}s). Last error: {e}"
                )


def fetch_seq_coordinates(genome: str, gggenome_url: str, seq_subset: str) -> dict:
    url = f"{gggenome_url}/{genome}/{seq_subset}"
    beds = fetch_bed_without_verification(url + ".bed")[1:]
    if len(beds) > 1:
        error_message = (
            f"The sequence matched multiple regions in the reference genome:\n"
            f"{gggenome_url}/{genome}/{seq_subset}\n\n"
            f"To resolve this, please specify the exact genomic coordinates by providing "
            f"a BED6 format file with the -b/--bed option.\n"
            f"See the documentation for details:\n"
            f"https://github.com/akikuno/DAJIN2#using-bed-files-for-genomic-coordinates"
        )
        raise ValueError(error_message)

    bed = beds[0]

    if bed == "### No items found. ###":
        error_message = (
            f"The sequence did not match any region in the reference genome:\n"
            f"{gggenome_url}/{genome}/{seq_subset}\n\n"
            f"To resolve this, please specify the exact genomic coordinates by providing "
            f"a BED6 format file with the -b/--bed option.\n"
            f"See the documentation for details:\n"
            f"https://github.com/akikuno/DAJIN2#using-bed-files-for-genomic-coordinates"
        )
        raise ValueError(error_message)

    chrom, start, end, *_, strand = bed

    return {"chrom": chrom, "strand": strand, "start": int(start), "end": int(end)}


def fetch_coordinates(genome: str, gggenome_url: str, seq: str) -> dict:
    seq_start, seq_end = seq[:100], seq[-100:]
    coordinate_start = fetch_seq_coordinates(genome, gggenome_url, seq_start)
    coordinate_end = fetch_seq_coordinates(genome, gggenome_url, seq_end)

    chromosome, strand = coordinate_start["chrom"], coordinate_start["strand"]
    if strand == "+":
        start, end = coordinate_start["start"], coordinate_end["end"]
    else:
        start, end = coordinate_end["start"], coordinate_start["end"]

    return {"genome": genome, "chrom": chromosome, "start": start, "end": end, "strand": strand}


def fetch_chromosome_size(genome: str, chrom: str, goldenpath_url: str) -> int:
    url = f"{goldenpath_url}/{genome}/bigZips/{genome}.chrom.sizes"

    records = fetch_bed_without_verification(url)
    for record in records:
        chrom_name, size = record
        if chrom == chrom_name:
            return int(size)

    raise ValueError(f"Chromosome {chrom} size not found.")
