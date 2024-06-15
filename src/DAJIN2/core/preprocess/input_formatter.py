from __future__ import annotations

import uuid
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from DAJIN2.core import preprocess
from DAJIN2.utils import config, fastx_handler, io


def parse_arguments(arguments: dict) -> tuple:
    genome_urls = defaultdict(str)
    if arguments.get("genome"):
        genome_urls.update(
            {"genome": arguments["genome"], "blat": arguments["blat"], "goldenpath": arguments["goldenpath"]}
        )

    return (
        Path(arguments["sample"]),
        Path(arguments["control"]),
        Path(arguments["allele"]),
        arguments["name"],
        arguments["threads"],
        genome_urls,
        uuid.uuid4().hex,
    )


def convert_input_paths_to_posix(sample: str, control: str, allele: str) -> tuple:
    sample = io.convert_to_posix(sample)
    control = io.convert_to_posix(control)
    allele = io.convert_to_posix(allele)

    return sample, control, allele


def create_temporal_directory(name: str, control_name: str) -> Path:
    tempdir = Path(config.TEMP_ROOT_DIR, name)
    Path(tempdir, "cache", ".igvjs", control_name).mkdir(parents=True, exist_ok=True)

    return tempdir


def check_caches(tempdir: Path, path_allele: str, genome_url: str) -> bool:
    is_cache_hash = preprocess.exists_cached_hash(tempdir=tempdir, path=path_allele)
    is_cache_genome = preprocess.exists_cached_genome(tempdir=tempdir, genome=genome_url)

    return is_cache_hash and is_cache_genome


def get_genome_coordinates(genome_urls: dict, fasta_alleles: dict, is_cache_genome: bool, tempdir: Path) -> dict:
    genome_coordinates = {
        "genome": genome_urls["genome"],
        "chrom_size": 0,
        "chrom": "control",
        "start": 0,
        "end": len(fasta_alleles["control"]) - 1,
        "strand": "+",
    }
    if genome_urls["genome"]:
        if is_cache_genome:
            genome_coordinates = next(io.read_jsonl(Path(tempdir, "cache", "genome_coordinates.jsonl")))
        else:
            genome_coordinates = preprocess.fetch_coordinates(
                genome_coordinates, genome_urls, fasta_alleles["control"]
            )
            genome_coordinates["chrom_size"] = preprocess.fetch_chromosome_size(genome_coordinates, genome_urls)
            io.write_jsonl([genome_coordinates], Path(tempdir, "cache", "genome_coordinates.jsonl"))

    return genome_coordinates


@dataclass(frozen=True)
class FormattedInputs:
    path_sample: Path
    path_control: Path
    path_allele: Path
    sample_name: str
    control_name: str
    fasta_alleles: dict[str, str]
    tempdir: Path
    genome_coordinates: dict[str, str]
    threads: int
    uuid: str


def format_inputs(arguments: dict) -> FormattedInputs:
    path_sample, path_control, path_allele, name, threads, genome_urls, uuid = parse_arguments(arguments)
    path_sample, path_control, path_allele = convert_input_paths_to_posix(path_sample, path_control, path_allele)
    sample_name = fastx_handler.extract_filename(path_sample)
    control_name = fastx_handler.extract_filename(path_control)
    fasta_alleles = fastx_handler.dictionize_allele(path_allele)
    tempdir = create_temporal_directory(name, control_name)
    is_cache_genome = check_caches(tempdir, path_allele, genome_urls["genome"])
    genome_coordinates = get_genome_coordinates(genome_urls, fasta_alleles, is_cache_genome, tempdir)

    return FormattedInputs(
        path_sample,
        path_control,
        path_allele,
        io.sanitize_name(sample_name),
        io.sanitize_name(control_name),
        fasta_alleles,
        tempdir,
        genome_coordinates,
        threads,
        uuid,
    )
