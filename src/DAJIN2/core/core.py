from __future__ import annotations

import os
import sys
import pickle

# import psutil
import resource
import shutil

from pathlib import Path
from datetime import datetime
from collections import defaultdict
from typing import NamedTuple

from DAJIN2.utils import io
from DAJIN2.core import classification, clustering, consensus, preprocess, report
from DAJIN2.utils.config import TEMP_ROOT_DIR

# # limit max memory usage
# # available_memory = psutil.virtual_memory().available
# available_memory = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
# _, hard_limit = resource.getrlimit(resource.RLIMIT_DATA)
# # if hard_limit == resource.RLIM_INFINITY:
# #     hard_limit = available_memory
# # Ensure that the new limit does not exceed the hard limit
# limit = min(int(available_memory * 0.9), hard_limit)
# resource.setrlimit(resource.RLIMIT_DATA, (limit, limit))


def parse_arguments(arguments: dict) -> tuple:
    genome_urls = defaultdict(str)
    if "genome" in arguments:
        genome_urls.update(
            {"genome": arguments["genome"], "blat": arguments["blat"], "goldenpath": arguments["goldenpath"]}
        )

    return (
        arguments["sample"],
        arguments["control"],
        arguments["allele"],
        arguments["name"],
        arguments["threads"],
        genome_urls,
    )


def convert_inputs_to_posix(sample: str, control: str, allele: str) -> tuple:
    sample = io.convert_to_posix(sample)
    control = io.convert_to_posix(control)
    allele = io.convert_to_posix(allele)
    return sample, control, allele


def create_temporal_directory(name: str, control_name: str) -> Path:
    tempdir = Path(TEMP_ROOT_DIR, name)
    Path(tempdir, "cache", ".igvjs", control_name).mkdir(parents=True, exist_ok=True)
    return tempdir


def check_caches(tempdir: Path, path_allele: str, genome_url: str) -> bool:
    is_cache_hash = preprocess.cache_checker.exists_cached_hash(tempdir=tempdir, path=path_allele)
    is_cache_genome = preprocess.cache_checker.exists_cached_genome(tempdir=tempdir, genome=genome_url)
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
            genome_coordinates = list(io.read_jsonl(Path(tempdir, "cache", "genome_coodinates.jsonl")))[0]
        else:
            genome_coordinates = preprocess.genome_fetcher.fetch_coordinates(
                genome_coordinates, genome_urls, fasta_alleles["control"]
            )
            genome_coordinates["chrom_size"] = preprocess.genome_fetcher.fetch_chromosome_size(
                genome_coordinates, genome_urls
            )
            io.write_jsonl([genome_coordinates], Path(tempdir, "cache", "genome_coodinates.jsonl"))
    return genome_coordinates


class FormattedInputs(NamedTuple):
    path_sample: str
    path_control: str
    path_allele: str
    sample_name: str
    control_name: str
    fasta_alleles: dict[str, str]
    tempdir: Path
    genome_coordinates: dict[str, str]
    threads: int


def format_inputs(arguments: dict) -> FormattedInputs:
    path_sample, path_control, path_allele, name, threads, genome_urls = parse_arguments(arguments)
    path_sample, path_control, path_allele = convert_inputs_to_posix(path_sample, path_control, path_allele)
    sample_name = preprocess.fastx_parser.extract_basename(path_sample)
    control_name = preprocess.fastx_parser.extract_basename(path_control)
    fasta_alleles = preprocess.fastx_parser.dictionize_allele(path_allele)
    tempdir = create_temporal_directory(name, control_name)
    is_cache_genome = check_caches(tempdir, path_allele, genome_urls["genome"])
    genome_coordinates = get_genome_coordinates(genome_urls, fasta_alleles, is_cache_genome, tempdir)
    return FormattedInputs(
        path_sample,
        path_control,
        path_allele,
        sample_name,
        control_name,
        fasta_alleles,
        tempdir,
        genome_coordinates,
        threads,
    )


def _dtnow() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


###########################################################
# main
###########################################################


def execute_control(arguments: dict):
    print(f"{_dtnow()}: {arguments['control']} is now processing...", file=sys.stderr)
    ###########################################################
    # Preprocess
    ###########################################################
    ARGS = format_inputs(arguments)
    preprocess.directories.create_temporal(ARGS.tempdir, ARGS.control_name, is_control=True)
    preprocess.directories.create_report(ARGS.tempdir, ARGS.control_name, is_control=True)
    io.cache_control_hash(ARGS.tempdir, ARGS.path_allele)
    ###########################################################
    # Check caches
    ###########################################################
    if Path(ARGS.tempdir, "report", "BAM", ARGS.control_name, f"{ARGS.control_name}.bam").exists():
        print(
            f"{arguments['control']} is already preprocessed and reuse the results for the current run...",
            file=sys.stderr,
        )
        return
    print(f"{_dtnow()}: Preprocess {arguments['control']}...", file=sys.stderr)
    ###########################################################
    # Mapping
    ###########################################################
    # ============================================================
    # Export fasta files as single-FASTA format
    # ============================================================
    preprocess.fastx_parser.export_fasta_files(ARGS.tempdir, ARGS.fasta_alleles, ARGS.control_name)
    # ============================================================
    # Mapping using mappy
    # ============================================================
    paths_fasta = Path(ARGS.tempdir, ARGS.control_name, "fasta").glob("*.fasta")
    preprocess.mapping.generate_sam(ARGS.tempdir, paths_fasta, ARGS.path_control, ARGS.control_name, ARGS.threads)
    ###########################################################
    # MIDSV conversion
    ###########################################################
    preprocess.midsv_caller.execute(ARGS.tempdir, ARGS.fasta_alleles, ARGS.control_name)
    ###########################################################
    # Prepare data to `extract mutaion loci`
    ###########################################################
    preprocess.extract_mutation_loci(
        ARGS.tempdir, ARGS.fasta_alleles, ARGS.sample_name, ARGS.control_name, is_control=True
    )
    ###########################################################
    # Output BAM
    ###########################################################
    print(f"{_dtnow()}: Output BAM files of {arguments['control']}...", file=sys.stderr)
    report.report_bam.output_bam(
        ARGS.tempdir, ARGS.control_name, ARGS.genome_coordinates, ARGS.threads, is_control=True
    )
    ###########################################################
    # Finish call
    ###########################################################
    print(f"{_dtnow()}: \N{teacup without handle} {arguments['control']} is finished!", file=sys.stderr)


def execute_sample(arguments: dict):
    print(f"{_dtnow()}: {arguments['sample']} is now processing...", file=sys.stderr)
    ###########################################################
    # Preprocess
    ###########################################################
    ARGS = format_inputs(arguments)
    preprocess.directories.create_temporal(ARGS.tempdir, ARGS.sample_name)
    preprocess.directories.create_report(ARGS.tempdir, ARGS.sample_name)

    print(f"{_dtnow()}: Preprocess {arguments['sample']}...", file=sys.stderr)

    for path_fasta in Path(ARGS.tempdir, ARGS.control_name, "fasta").glob("*.fasta"):
        shutil.copy(path_fasta, Path(ARGS.tempdir, ARGS.sample_name, "fasta"))
    # ============================================================
    # Mapping with mappy
    # ============================================================
    paths_fasta = Path(ARGS.tempdir, ARGS.sample_name, "fasta").glob("*.fasta")
    preprocess.mapping.generate_sam(ARGS.tempdir, paths_fasta, ARGS.path_sample, ARGS.sample_name, ARGS.threads)
    # ============================================================
    # MIDSV conversion
    # ============================================================
    preprocess.midsv_caller.execute(ARGS.tempdir, ARGS.fasta_alleles, ARGS.sample_name)
    # ============================================================
    # Extract mutation loci
    # ============================================================
    preprocess.extract_knockin_loci(ARGS.tempdir, ARGS.sample_name)
    preprocess.extract_mutation_loci(ARGS.tempdir, ARGS.fasta_alleles, ARGS.sample_name, ARGS.control_name)
    # ============================================================
    # Detect and align insertion alleles
    # ============================================================
    paths_predifined_allele = {str(p) for p in Path(ARGS.tempdir, ARGS.sample_name, "fasta").glob("*.fasta")}
    preprocess.generate_insertion_fasta(ARGS.tempdir, ARGS.sample_name, ARGS.control_name, ARGS.fasta_alleles)
    paths_insertion = {str(p) for p in Path(ARGS.tempdir, ARGS.sample_name, "fasta").glob("insertion*.fasta")}
    paths_insertion -= paths_predifined_allele
    if paths_insertion:
        # mapping to insertion alleles
        preprocess.mapping.generate_sam(
            ARGS.tempdir, paths_insertion, ARGS.path_control, ARGS.control_name, ARGS.threads
        )
        preprocess.mapping.generate_sam(ARGS.tempdir, paths_insertion, ARGS.path_sample, ARGS.sample_name, ARGS.threads)
        # add insertions to ARGS.fasta_alleles
        for path_fasta in paths_insertion:
            allele, seq = Path(path_fasta).read_text().strip().split("\n")
            allele = allele.replace(">", "")
            ARGS.fasta_alleles[allele] = seq
        # MIDSV conversion
        preprocess.midsv_caller.execute(ARGS.tempdir, ARGS.fasta_alleles, ARGS.control_name)
        preprocess.midsv_caller.execute(ARGS.tempdir, ARGS.fasta_alleles, ARGS.sample_name)
        # Reculculate mutation loci
        preprocess.extract_mutation_loci(
            ARGS.tempdir, ARGS.fasta_alleles, ARGS.sample_name, ARGS.control_name, is_control=True
        )
        preprocess.extract_knockin_loci(ARGS.tempdir, ARGS.sample_name)
        preprocess.extract_mutation_loci(ARGS.tempdir, ARGS.fasta_alleles, ARGS.sample_name, ARGS.control_name)
    ########################################################################
    # Classify alleles
    ########################################################################
    print(f"{_dtnow()}: Classify {arguments['sample']}...", file=sys.stderr)
    classif_sample = classification.classify_alleles(ARGS.tempdir, ARGS.fasta_alleles, ARGS.sample_name)
    with open(Path(ARGS.tempdir, ARGS.sample_name, "classif_sample.pickle"), "wb") as p:
        pickle.dump(classif_sample, p)
    ########################################################################
    # Clustering
    ########################################################################
    print(f"{_dtnow()}: Clustering {arguments['sample']}...", file=sys.stderr)
    clust_sample = clustering.add_labels(classif_sample, ARGS.tempdir, ARGS.sample_name, ARGS.control_name)
    clust_sample = clustering.add_readnum(clust_sample)
    clust_sample = clustering.add_percent(clust_sample)
    clust_sample = clustering.update_labels(clust_sample)
    with open(Path(ARGS.tempdir, ARGS.sample_name, "clust_sample.pickle"), "wb") as p:
        pickle.dump(clust_sample, p)
    ########################################################################
    # Consensus call
    ########################################################################
    print(f"{_dtnow()}: Consensus calling of {arguments['sample']}...", file=sys.stderr)
    # Downsampling to 1000 reads in each LABEL
    clust_subset_sample = consensus.subset_clust(clust_sample, 1000)
    cons_percentage, cons_sequence = consensus.call_consensus(ARGS.tempdir, ARGS.sample_name, clust_subset_sample)
    # cons_percentage, cons_sequence = consensus.call_consensus(clust_subset_sample, MUTATION_LOCI_LABELS)
    allele_names = consensus.call_allele_name(cons_sequence, cons_percentage, ARGS.fasta_alleles)
    cons_percentage = consensus.update_key_by_allele_name(cons_percentage, allele_names)
    cons_sequence = consensus.update_key_by_allele_name(cons_sequence, allele_names)
    RESULT_SAMPLE = consensus.add_key_by_allele_name(clust_sample, allele_names)
    RESULT_SAMPLE.sort(key=lambda x: x["LABEL"])
    ########################################################################
    # Output Report：RESULT/FASTA/HTML/BAM
    ########################################################################
    print(f"{_dtnow()}: Output reports of {arguments['sample']}...", file=sys.stderr)
    # RESULT
    io.write_jsonl(RESULT_SAMPLE, Path(ARGS.tempdir, "result", f"{ARGS.sample_name}.jsonl"))
    # FASTA
    report.report_files.to_fasta(ARGS.tempdir, ARGS.sample_name, cons_sequence)
    # HTML
    report.report_files.to_html(ARGS.tempdir, ARGS.sample_name, cons_percentage)
    # CSV (Allele Info)
    report.report_mutation.to_csv(ARGS.tempdir, ARGS.sample_name, ARGS.genome_coordinates, cons_percentage)
    # BAM
    report.report_bam.output_bam(ARGS.tempdir, ARGS.sample_name, ARGS.genome_coordinates, ARGS.threads, RESULT_SAMPLE)
    for path_bam_igvjs in Path(ARGS.tempdir, "cache", ".igvjs").glob(f"{ARGS.control_name}_control.bam*"):
        shutil.copy(path_bam_igvjs, Path(ARGS.tempdir, "report", ".igvjs", ARGS.sample_name))
    # VCF
    # working in progress
    ###########################################################
    # Finish call
    ###########################################################
    print(f"{_dtnow()}: \N{teacup without handle} {arguments['sample']} is finished!", file=sys.stderr)