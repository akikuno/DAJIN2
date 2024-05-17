from __future__ import annotations

import ssl
from urllib.request import urlopen


def fetch_html_without_verification(url: str) -> str:
    context = ssl._create_unverified_context()  # Create an SSL context that temporarily disables verification
    with urlopen(url, context=context, timeout=10) as response:
        return response.read().decode("utf-8").split("\n")


def fetch_seq_coordinates(genome: str, blat_url: str, seq: str) -> dict:
    url = f"{blat_url}?db={genome}&type=BLAT&userSeq={seq}"
    records = fetch_html_without_verification(url)
    matches = []
    for record in records:
        if "100.0%" not in record:
            continue
        record_trim = [r for r in record.split(" ") if r]
        if record_trim[-1] == str(len(seq)):
            matches = record_trim

    if not matches:
        raise ValueError(f"{seq[:60]}... is not found in {genome}")

    chrom, strand, start, end, _ = matches[-5:]
    return {"chrom": chrom, "strand": strand, "start": int(start), "end": int(end)}


def fetch_coordinates(genome_coordinates: dict, genome_urls: dict, seq: str) -> dict:
    genome = genome_coordinates["genome"]
    blat_url = genome_urls["blat"]

    seq_start, seq_end = seq[:1000], seq[-1000:]
    coordinate_start = fetch_seq_coordinates(genome, blat_url, seq_start)
    coordinate_end = fetch_seq_coordinates(genome, blat_url, seq_end)

    chromosome, strand = coordinate_start["chrom"], coordinate_start["strand"]
    if strand == "+":
        start, end = coordinate_start["start"], coordinate_end["end"]
    else:
        start, end = coordinate_end["start"], coordinate_start["end"]

    return {"genome": genome, "chrom": chromosome, "start": start, "end": end, "strand": strand}


def fetch_chromosome_size(genome_coordinates: dict, genome_urls: dict) -> int:
    chrom = genome_coordinates["chrom"]
    genome = genome_coordinates["genome"]
    url = f"{genome_urls['goldenpath']}/{genome}/bigZips/{genome}.chrom.sizes"

    records = fetch_html_without_verification(url)
    for record in records:
        chrom_name, size = record.split("\t")
        if chrom == chrom_name:
            return int(size)
    raise ValueError(f"Chromosome {chrom} size not found.")
