from __future__ import annotations

import os
import sys
import pickle
import resource
import shutil
from datetime import datetime
from pathlib import Path

import midsv

from collections import defaultdict
from DAJIN2.utils import io
from DAJIN2.core import classification, clustering, consensus, preprocess, report
from DAJIN2.utils.config import TEMP_ROOT_DIR

# limit max memory usage
mem_bytes = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
resource.setrlimit(resource.RLIMIT_DATA, (int(mem_bytes * 9 / 10), -1))


def parse_arguments(arguments: dict):
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


def check_caches(control: str, tempdir: Path, genome_url: str) -> tuple:
    is_cache_control = preprocess.check_caches.exists_cached_control(control, tempdir)
    is_cache_genome = preprocess.check_caches.exists_cached_genome(genome_url, tempdir, is_cache_control)
    return is_cache_control, is_cache_genome


def get_genome_coordinates(genome_urls: dict, fasta_alleles: dict, is_cache_genome: bool, tempdir: Path) -> dict:
    genome_coordinates = {
        "genome": genome_urls["genome"],
        "chrom_size": 0,
        "chr": "control",
        "start": 0,
        "end": len(fasta_alleles["control"]) - 1,
        "strand": "+",
    }
    if genome_urls["genome"] and not is_cache_genome:
        genome_coordinates = preprocess.format_inputs.fetch_coordinate(
            genome_coordinates, genome_urls, fasta_alleles["control"]
        )
        genome_coordinates = preprocess.format_inputs.fetch_chrom_size(genome_coordinates, genome_urls)
        midsv.write_jsonl([genome_coordinates], Path(tempdir, "cache", "genome_coodinates.jsonl"))
    elif genome_urls["genome"]:
        genome_coordinates = midsv.read_jsonl(Path(tempdir, "cache", "genome_coodinates.jsonl"))
    return genome_coordinates


def format_inputs(arguments: dict) -> tuple:
    sample, control, allele, name, threads, genome_urls = parse_arguments(arguments)
    sample, control, allele = convert_inputs_to_posix(sample, control, allele)
    sample_name = preprocess.format_inputs.extract_basename(sample)
    control_name = preprocess.format_inputs.extract_basename(control)
    fasta_alleles = preprocess.format_inputs.dictionize_allele(allele)
    tempdir = create_temporal_directory(name, control_name)
    is_cache_control, is_cache_genome = check_caches(control, tempdir, genome_urls["genome"])
    genome_coordinates = get_genome_coordinates(genome_urls, fasta_alleles, is_cache_genome, tempdir)
    return sample_name, control_name, fasta_alleles, tempdir, genome_coordinates, threads


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
    SAMPLE, CONTROL, ALLELE, NAME, THREADS, GENOME_URLS = parse_arguments(arguments)
    SAMPLE_NAME, CONTROL_NAME, FASTA_ALLELES, TEMPDIR, GENOME_COODINATES, THREADS = format_inputs(arguments)
    preprocess.format_inputs.make_directories(TEMPDIR, CONTROL_NAME, is_control=True)
    preprocess.format_inputs.make_report_directories(TEMPDIR, CONTROL_NAME, is_control=True)
    ###########################################################
    # Check caches
    ###########################################################
    if Path(TEMPDIR, "report", "BAM", CONTROL_NAME, f"{CONTROL_NAME}.bam").exists():
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
    preprocess.format_inputs.export_fasta_files(TEMPDIR, FASTA_ALLELES, CONTROL_NAME)
    # ============================================================
    # Mapping using mappy
    # ============================================================
    paths_fasta = Path(TEMPDIR, CONTROL_NAME, "fasta").glob("*.fasta")
    preprocess.align.generate_sam(TEMPDIR, paths_fasta, CONTROL, CONTROL_NAME, THREADS)
    ###########################################################
    # MIDSV conversion
    ###########################################################
    preprocess.call_midsv(TEMPDIR, FASTA_ALLELES, CONTROL_NAME)
    ###########################################################
    # Prepare data to `extract mutaion loci`
    ###########################################################
    # preprocess.save_index_mapping(TEMPDIR)
    preprocess.extract_mutation_loci(TEMPDIR, FASTA_ALLELES, SAMPLE_NAME, CONTROL_NAME, is_control=True)
    ###########################################################
    # Output BAM
    ###########################################################
    print(f"{_dtnow()}: Output BAM files of {arguments['control']}...", file=sys.stderr)
    report.report_bam.output_bam(TEMPDIR, CONTROL_NAME, GENOME_COODINATES, THREADS, is_control=True)
    ###########################################################
    # Finish call
    ###########################################################
    print(f"{_dtnow()}: \N{teacup without handle} {arguments['control']} is finished!", file=sys.stderr)


def execute_sample(arguments: dict):
    print(f"{_dtnow()}: {arguments['sample']} is now processing...", file=sys.stderr)
    ###########################################################
    # Preprocess
    ###########################################################
    SAMPLE, CONTROL, ALLELE, NAME, THREADS, GENOME_URLS = parse_arguments(arguments)
    SAMPLE_NAME, CONTROL_NAME, FASTA_ALLELES, TEMPDIR, GENOME_COODINATES, THREADS = format_inputs(arguments)
    preprocess.format_inputs.make_directories(TEMPDIR, SAMPLE_NAME)
    preprocess.format_inputs.make_report_directories(TEMPDIR, SAMPLE_NAME)

    print(f"{_dtnow()}: Preprocess {arguments['sample']}...", file=sys.stderr)

    for path_fasta in Path(TEMPDIR, CONTROL_NAME, "fasta").glob("*.fasta"):
        shutil.copy(path_fasta, Path(TEMPDIR, SAMPLE_NAME, "fasta"))
    # ============================================================
    # Mapping with mappy
    # ============================================================
    paths_fasta = Path(TEMPDIR, SAMPLE_NAME, "fasta").glob("*.fasta")
    preprocess.align.generate_sam(TEMPDIR, paths_fasta, SAMPLE, SAMPLE_NAME, THREADS)
    # ============================================================
    # MIDSV conversion
    # ============================================================
    preprocess.call_midsv(TEMPDIR, FASTA_ALLELES, SAMPLE_NAME)
    # ============================================================
    # Extract mutation loci
    # ============================================================
    preprocess.extract_knockin_loci(TEMPDIR, SAMPLE_NAME)
    preprocess.extract_mutation_loci(TEMPDIR, FASTA_ALLELES, SAMPLE_NAME, CONTROL_NAME)
    # ============================================================
    # Detect and align insertion alleles
    # ============================================================
    paths_predifined_allele = {str(p) for p in Path(TEMPDIR, SAMPLE_NAME, "fasta").glob("*.fasta")}
    preprocess.generate_insertion_fasta(TEMPDIR, SAMPLE_NAME, CONTROL_NAME, FASTA_ALLELES)
    paths_insertion = {str(p) for p in Path(TEMPDIR, SAMPLE_NAME, "fasta").glob("insertion*.fasta")}
    paths_insertion -= paths_predifined_allele
    if paths_insertion:
        # mapping to insertion alleles
        preprocess.align.generate_sam(TEMPDIR, paths_insertion, CONTROL, CONTROL_NAME, THREADS)
        preprocess.align.generate_sam(TEMPDIR, paths_insertion, SAMPLE, SAMPLE_NAME, THREADS)
        # add insertions to FASTA_ALLELES
        for path_fasta in paths_insertion:
            allele, seq = Path(path_fasta).read_text().strip().split("\n")
            allele = allele.replace(">", "")
            FASTA_ALLELES[allele] = seq
        # MIDSV conversion
        preprocess.call_midsv(TEMPDIR, FASTA_ALLELES, CONTROL_NAME)
        preprocess.call_midsv(TEMPDIR, FASTA_ALLELES, SAMPLE_NAME)
        # Reculculate mutation loci
        preprocess.extract_mutation_loci(TEMPDIR, FASTA_ALLELES, SAMPLE_NAME, CONTROL_NAME, is_control=True)
        preprocess.extract_knockin_loci(TEMPDIR, SAMPLE_NAME)
        preprocess.extract_mutation_loci(TEMPDIR, FASTA_ALLELES, SAMPLE_NAME, CONTROL_NAME)
    ########################################################################
    # Classify alleles
    ########################################################################
    print(f"{_dtnow()}: Classify {arguments['sample']}...", file=sys.stderr)
    classif_sample = classification.classify_alleles(TEMPDIR, FASTA_ALLELES, SAMPLE_NAME)
    with open(Path(TEMPDIR, SAMPLE_NAME, "classif_sample.pickle"), "wb") as p:
        pickle.dump(classif_sample, p)
    ########################################################################
    # Clustering
    ########################################################################
    print(f"{_dtnow()}: Clustering {arguments['sample']}...", file=sys.stderr)
    clust_sample = clustering.add_labels(classif_sample, TEMPDIR, SAMPLE_NAME, CONTROL_NAME, THREADS)
    clust_sample = clustering.add_readnum(clust_sample)
    clust_sample = clustering.add_percent(clust_sample)
    clust_sample = clustering.update_labels(clust_sample)
    with open(Path(TEMPDIR, SAMPLE_NAME, "clust_sample.pickle"), "wb") as p:
        pickle.dump(clust_sample, p)
    ########################################################################
    # Consensus call
    ########################################################################
    print(f"{_dtnow()}: Consensus calling of {arguments['sample']}...", file=sys.stderr)
    # Downsampling to 1000 reads in each LABEL
    clust_subset_sample = consensus.subset_clust(clust_sample, 1000)
    cons_percentage, cons_sequence = consensus.call_consensus(TEMPDIR, SAMPLE_NAME, clust_subset_sample)
    # cons_percentage, cons_sequence = consensus.call_consensus(clust_subset_sample, MUTATION_LOCI_LABELS)
    allele_names = consensus.call_allele_name(cons_sequence, cons_percentage, FASTA_ALLELES)
    cons_percentage = consensus.update_key_by_allele_name(cons_percentage, allele_names)
    cons_sequence = consensus.update_key_by_allele_name(cons_sequence, allele_names)
    RESULT_SAMPLE = consensus.add_key_by_allele_name(clust_sample, allele_names)
    RESULT_SAMPLE.sort(key=lambda x: x["LABEL"])
    ########################################################################
    # Output Report：RESULT/FASTA/HTML/BAM
    ########################################################################
    print(f"{_dtnow()}: Output reports of {arguments['sample']}...", file=sys.stderr)
    # RESULT
    midsv.write_jsonl(RESULT_SAMPLE, Path(TEMPDIR, "result", f"{SAMPLE_NAME}.jsonl"))
    # FASTA
    report.report_files.to_fasta(TEMPDIR, SAMPLE_NAME, cons_sequence)
    # HTML
    report.report_files.to_html(TEMPDIR, SAMPLE_NAME, cons_percentage)
    # CSV (Allele Info)
    report.report_mutation.to_csv(TEMPDIR, SAMPLE_NAME, GENOME_COODINATES, cons_percentage)
    # BAM
    report.report_bam.output_bam(TEMPDIR, SAMPLE_NAME, GENOME_COODINATES, THREADS, RESULT_SAMPLE)
    for path_bam_igvjs in Path(TEMPDIR, "cache", ".igvjs").glob(f"{CONTROL_NAME}_control.bam*"):
        shutil.copy(path_bam_igvjs, Path(TEMPDIR, "report", ".igvjs", SAMPLE_NAME))
    # VCF
    # working in progress
    ###########################################################
    # Finish call
    ###########################################################
    print(f"{_dtnow()}: \N{teacup without handle} {arguments['sample']} is finished!", file=sys.stderr)
