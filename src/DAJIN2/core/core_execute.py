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
from DAJIN2.core import classification, clustering, consensus, preprocess, report

# limit max memory usage
mem_bytes = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
resource.setrlimit(resource.RLIMIT_DATA, (int(mem_bytes * 9 / 10), -1))


def _parse_arguments(arguments: dict):
    SAMPLE: str = arguments["sample"]
    CONTROL: str = arguments["control"]
    ALLELE: str = arguments["allele"]
    NAME: str = arguments["name"]
    THREADS: int = arguments["threads"]
    GENOME_URLS = defaultdict(str)
    if "genome" in arguments:
        GENOME_URLS["genome"] = arguments["genome"]
        GENOME_URLS["blat"] = arguments["blat"]
        GENOME_URLS["goldenpath"] = arguments["goldenpath"]
    return SAMPLE, CONTROL, ALLELE, NAME, THREADS, GENOME_URLS


def _format_inputs(arguments: dict):
    SAMPLE, CONTROL, ALLELE, NAME, THREADS, GENOME_URLS = _parse_arguments(arguments)
    SAMPLE = preprocess.format_inputs.convert_to_posix_path(SAMPLE)
    CONTROL = preprocess.format_inputs.convert_to_posix_path(CONTROL)
    ALLELE = preprocess.format_inputs.convert_to_posix_path(ALLELE)

    SAMPLE_NAME: str = preprocess.format_inputs.extract_basename(SAMPLE)
    CONTROL_NAME: str = preprocess.format_inputs.extract_basename(CONTROL)
    FASTA_ALLELES: dict = preprocess.format_inputs.dictionize_allele(ALLELE)

    TEMPDIR = Path("DAJINResults", ".tempdir", NAME)
    Path(TEMPDIR, "cache", ".igvjs", CONTROL_NAME).mkdir(parents=True, exist_ok=True)

    IS_CACHE_CONTROL = preprocess.check_caches.exists_cached_control(CONTROL, TEMPDIR)
    IS_CACHE_GENOME = preprocess.check_caches.exists_cached_genome(GENOME_URLS["genome"], TEMPDIR, IS_CACHE_CONTROL)

    GENOME_COODINATES = {
        "genome": GENOME_URLS["genome"],
        "chrom_size": 0,
        "chr": "control",
        "start": 0,
        "end": len(FASTA_ALLELES["control"]) - 1,
        "strand": "+",
    }
    if GENOME_URLS["genome"]:
        if not IS_CACHE_GENOME:
            GENOME_COODINATES = preprocess.format_inputs.fetch_coordinate(
                GENOME_COODINATES, GENOME_URLS, FASTA_ALLELES["control"]
            )
            GENOME_COODINATES = preprocess.format_inputs.fetch_chrom_size(GENOME_COODINATES, GENOME_URLS)
            midsv.write_jsonl([GENOME_COODINATES], Path(TEMPDIR, "cache", "genome_coodinates.jsonl"))
        else:
            GENOME_COODINATES = midsv.read_jsonl(Path(TEMPDIR, "cache", "genome_coodinates.jsonl"))
    return SAMPLE_NAME, CONTROL_NAME, FASTA_ALLELES, TEMPDIR, GENOME_COODINATES, THREADS


def _dtnow() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def execute_control(arguments: dict):
    print(f"{_dtnow()}: {arguments['control']} is now processing...", file=sys.stderr)
    ###########################################################
    # Preprocess
    ###########################################################
    SAMPLE, CONTROL, ALLELE, NAME, THREADS, GENOME_URLS = _parse_arguments(arguments)
    SAMPLE_NAME, CONTROL_NAME, FASTA_ALLELES, TEMPDIR, GENOME_COODINATES, THREADS = _format_inputs(arguments)
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
    SAMPLE, CONTROL, ALLELE, NAME, THREADS, GENOME_URLS = _parse_arguments(arguments)
    SAMPLE_NAME, CONTROL_NAME, FASTA_ALLELES, TEMPDIR, GENOME_COODINATES, THREADS = _format_inputs(arguments)
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
    # MUTATION_LOCI_LABELS = consensus.extract_mutation_loci_by_labels(
    #     clust_sample, TEMPDIR, FASTA_ALLELES, CONTROL_NAME, SAMPLE_NAME
    # )
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
