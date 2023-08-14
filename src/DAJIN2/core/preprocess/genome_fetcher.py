from __future__ import annotations

from urllib.request import urlopen


def _fetch_seq_coordinates(genome: str, blat_url: str, seq: str) -> dict:
    url = f"{blat_url}?db={genome}&type=BLAT&userSeq={seq}"
    response = urlopen(url).read().decode("utf8").split("\n")
    matches = [x for x in response if "100.0%" in x]
    if not matches:
        raise ValueError(f"{seq} is not found in {genome}")
    chrom, strand, start, end, _ = matches[0].split()[-5:]
    return {"chrom": chrom, "strand": strand, "start": int(start), "end": int(end)}


def fetch_coordinates(genome_coordinates: dict, genome_urls: dict, seq: str) -> dict:
    genome = genome_coordinates["genome"]
    blat_url = genome_urls["blat"]

    seq_start, seq_end = seq[:1000], seq[-1000:]
    coordinate_start = _fetch_seq_coordinates(genome, blat_url, seq_start)
    coordinate_end = _fetch_seq_coordinates(genome, blat_url, seq_end)

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

    response = urlopen(url).read().decode("utf8").split("\n")
    for line in response:
        chrom_name, size = line.split("\t")
        if chrom == chrom_name:
            return int(size)
    raise ValueError(f"Chromosome {chrom} size not found.")
