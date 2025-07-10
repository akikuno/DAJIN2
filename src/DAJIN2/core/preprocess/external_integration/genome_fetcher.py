from __future__ import annotations

import re
import ssl
from urllib.request import urlopen


def fetch_html_without_verification(url: str) -> list[str]:
    context = ssl._create_unverified_context()  # Create an SSL context that temporarily disables verification
    with urlopen(url, context=context, timeout=10) as response:
        return response.read().decode("utf-8").split("\n")


def extract_blat_result(records: list[str], genome: str) -> dict:
    for record in records:
        if "YourSeq" in record and "details" in record:
            break
    record_trim = [r for r in record.split(" ") if r]

    chrom, strand, start, end, _ = record_trim[-5:]
    return {"chrom": chrom, "strand": strand, "start": int(start), "end": int(end), "genome": genome}


def extract_hub_id(records: list[str], genome: str) -> str:
    for record in records:
        if "db=hub_" in record:
            match = re.search(rf"db=(hub_\d+_{genome})", record)
            if match:
                return match.group(1)
    return ""


def fetch_seq_coordinates(genome: str, blat_url: str, seq_all: str, seq_subset: str, genome_urls: dict) -> dict:
    url = f"{blat_url}?db={genome}&type=BLAT&userSeq={seq_subset}"
    records = fetch_html_without_verification(url)
    matches = []
    for record in records:
        if "100.0%" not in record:
            continue
        record_trim = [r for r in record.split(" ") if r]
        if record_trim[-1] == str(len(seq_subset)):
            matches = record_trim

    if not matches:
        url = f"{blat_url}?db={genome}&type=BLAT&userSeq={seq_all}"
        records = fetch_html_without_verification(url)
        coords = extract_blat_result(records, genome)
        hub_id = extract_hub_id(records, genome)
        chrom_size = fetch_chromosome_size(coords, genome_urls)
        error_message = (
            f"There are differences with the UCSC reference sequence. "
            f"Please compare the DNA sequence in the following reference genome region with your FASTA file sequence to verify: "
            f"https://genome.ucsc.edu/cgi-bin/hgTracks?db={hub_id}&position={coords['chrom']}:{coords['start']}-{coords['end']}\n\n"
            f"Alternative: Use BED6 format file with -b/--bed option to specify coordinates directly.\n"
            f"{coords['chrom']}\t{coords['start']}\t{coords['end']}\t{chrom_size}\t0\t{coords['strand']}"
        )
        raise ValueError(error_message)

    chrom, strand, start, end, _ = matches[-5:]
    return {"chrom": chrom, "strand": strand, "start": int(start), "end": int(end)}


def fetch_coordinates(genome_coordinates: dict, genome_urls: dict, seq: str) -> dict:
    genome = genome_coordinates["genome"]
    blat_url = genome_urls["blat"]

    seq_start, seq_end = seq[:1000], seq[-1000:]
    coordinate_start = fetch_seq_coordinates(genome, blat_url, seq, seq_start, genome_urls)
    coordinate_end = fetch_seq_coordinates(genome, blat_url, seq, seq_end, genome_urls)

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
