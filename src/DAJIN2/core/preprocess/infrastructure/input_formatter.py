from __future__ import annotations

import json
import uuid
from dataclasses import dataclass
from pathlib import Path

from DAJIN2.core.preprocess.genome_coordinate.genome_coordinate_generator import (
    get_genome_coordinates_from_bed,
    get_genome_coordinates_from_server,
)
from DAJIN2.utils import config, fastx_handler, io


def parse_arguments(arguments: dict, tempdir: Path) -> tuple:
    genome = arguments.get("genome")
    bed = arguments.get("genome_coordinate")
    cache_file = Path(tempdir, "cache", "genome_coordinates.jsonl")
    cache_file.parent.mkdir(parents=True, exist_ok=True)

    genome_coordinates = None

    def fetch_and_cache_genome_coordinates(
        genome: str, gggenome_url: str, goldenpath_url: str, control_seq: str
    ) -> dict:
        gc = get_genome_coordinates_from_server(genome, gggenome_url, goldenpath_url, control_seq)
        with open(cache_file, "w") as f:
            json.dump(gc, f)
            f.write("\n")
        return gc

    if genome:
        gggenome_url = arguments["gggenome"]
        goldenpath_url = arguments["goldenpath"]
        control_seq = fastx_handler.dictionize_allele(arguments["allele"])["control"]
        if cache_file.exists():
            genome_coordinates = next(io.read_jsonl(cache_file))
            if genome_coordinates.get("genome") != genome:
                genome_coordinates = fetch_and_cache_genome_coordinates(
                    genome, gggenome_url, goldenpath_url, control_seq
                )
        else:
            genome_coordinates = fetch_and_cache_genome_coordinates(genome, gggenome_url, goldenpath_url, control_seq)
    elif bed:
        genome_coordinates = get_genome_coordinates_from_bed(bed)

    return (
        Path(arguments["sample"]),
        Path(arguments["control"]),
        Path(arguments["allele"]),
        arguments["threads"],
        genome_coordinates,
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
    no_filter: bool = False


def format_inputs(arguments: dict) -> FormattedInputs:
    tempdir = create_temporal_directory(arguments["name"], arguments["control"])
    path_sample, path_control, path_allele, threads, genome_coordinates, uuid = parse_arguments(arguments, tempdir)
    path_sample, path_control, path_allele = convert_input_paths_to_posix(path_sample, path_control, path_allele)
    sample_name = fastx_handler.extract_filename(path_sample)
    control_name = fastx_handler.extract_filename(path_control)
    fasta_alleles = fastx_handler.dictionize_allele(path_allele)

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
        arguments.get("no_filter", False),
    )
