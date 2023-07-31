from __future__ import annotations

import os
import sys
import pickle
import resource
import shutil
from datetime import datetime
from pathlib import Path

import midsv

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
    if "genome" in arguments:
        GENOME: str = arguments["genome"]
        GENOME_URLS = dict()
        GENOME_URLS["blat"] = arguments["blat"]
        GENOME_URLS["goldenpath"] = arguments["goldenpath"]
    else:
        GENOME = ""
        GENOME_URLS = dict()
        GENOME_URLS["blat"] = ""
        GENOME_URLS["goldenpath"] = ""
    return SAMPLE, CONTROL, ALLELE, NAME, THREADS, GENOME, GENOME_URLS


def _format_inputs(arguments: dict):
    SAMPLE, CONTROL, ALLELE, NAME, THREADS, GENOME, GENOME_URLS = _parse_arguments(arguments)
    SAMPLE = preprocess.format_inputs.convert_to_posix_path(SAMPLE)
    CONTROL = preprocess.format_inputs.convert_to_posix_path(CONTROL)
    ALLELE = preprocess.format_inputs.convert_to_posix_path(ALLELE)

    SAMPLE_NAME: str = preprocess.format_inputs.extract_basename(SAMPLE)
    CONTROL_NAME: str = preprocess.format_inputs.extract_basename(CONTROL)
    FASTA_ALLELES: dict = preprocess.format_inputs.dictionize_allele(ALLELE)

    TEMPDIR = Path("DAJINResults", ".tempdir", NAME)
    SUBDIRS = [
        "cache",
        "fasta",
        "sam",
        "midsv",
        "classification",
        "clustering",
        "report",
        "result",
        "mutation_loci",
        "knockin_loci",
    ]
    SUBDIRS_REPORT = ["HTML", "FASTA", "BAM", "MUTATION_INFO", ".igvjs"]
    preprocess.format_inputs.make_directories(TEMPDIR, SUBDIRS, SUBDIRS_REPORT, SAMPLE_NAME, CONTROL_NAME)

    IS_CACHE_CONTROL = preprocess.check_caches.exists_cached_control(CONTROL, TEMPDIR)
    IS_CACHE_GENOME = preprocess.check_caches.exists_cached_genome(GENOME, TEMPDIR, IS_CACHE_CONTROL)

    GENOME_COODINATES = {
        "genome": GENOME,
        "chrom_size": 0,
        "chr": "control",
        "start": 0,
        "end": len(FASTA_ALLELES["control"]) - 1,
        "strand": "+",
    }
    if GENOME:
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
    SAMPLE, CONTROL, ALLELE, NAME, THREADS, GENOME, GENOME_URLS = _parse_arguments(arguments)
    _, CONTROL_NAME, FASTA_ALLELES, TEMPDIR, GENOME_COODINATES, THREADS = _format_inputs(arguments)
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
    for identifier, sequence in FASTA_ALLELES.items():
        contents = "\n".join([">" + identifier, sequence]) + "\n"
        output_fasta = Path(TEMPDIR, "fasta", f"{identifier}.fasta")
        output_fasta.write_text(contents)
    # ============================================================
    # Mapping using mappy
    # ============================================================
    for path_fasta in Path(TEMPDIR, "fasta").glob("*.fasta"):
        name_fasta = path_fasta.stem
        preprocess.mappy_align.output_sam(TEMPDIR, path_fasta, name_fasta, CONTROL, CONTROL_NAME, threads=THREADS)
        preprocess.mappy_align.output_sam(
            TEMPDIR, path_fasta, name_fasta, CONTROL, CONTROL_NAME, preset="splice", threads=THREADS
        )
    ###########################################################
    # MIDSV conversion
    ###########################################################
    preprocess.call_midsv(TEMPDIR, FASTA_ALLELES, CONTROL_NAME)
    ###########################################################
    # Prepare data to `extract mutaion loci`
    ###########################################################
    preprocess.save_index_mapping(TEMPDIR)
    preprocess.process_mutation_loci(TEMPDIR, FASTA_ALLELES, CONTROL_NAME)
    ###########################################################
    # Output BAM
    ###########################################################
    print(f"{_dtnow()}: Output BAM files of {arguments['control']}...", file=sys.stderr)
    report.report_bam.output_bam_control(TEMPDIR, CONTROL_NAME, GENOME, GENOME_COODINATES, THREADS)
    ###########################################################
    # Finish call
    ###########################################################
    print(f"{_dtnow()}: \N{teacup without handle} {arguments['control']} is finished!", file=sys.stderr)


def execute_sample(arguments: dict):
    print(f"{_dtnow()}: {arguments['sample']} is now processing...", file=sys.stderr)
    ###########################################################
    # Preprocess
    ###########################################################
    SAMPLE, CONTROL, ALLELE, NAME, THREADS, GENOME, GENOME_URLS = _parse_arguments(arguments)
    SAMPLE_NAME, CONTROL_NAME, FASTA_ALLELES, TEMPDIR, GENOME_COODINATES, THREADS = _format_inputs(arguments)
    print(f"{_dtnow()}: Preprocess {arguments['sample']}...", file=sys.stderr)
    # ============================================================
    # Mapping with mappy
    # ============================================================
    for path_fasta in Path(TEMPDIR, "fasta").glob("*.fasta"):
        name_fasta = path_fasta.stem
        preprocess.mappy_align.output_sam(TEMPDIR, path_fasta, name_fasta, SAMPLE, SAMPLE_NAME, threads=THREADS)
        preprocess.mappy_align.output_sam(
            TEMPDIR, path_fasta, name_fasta, SAMPLE, SAMPLE_NAME, preset="splice", threads=THREADS
        )
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
    preprocess.generate_insertion_fasta(TEMPDIR, SAMPLE_NAME, CONTROL_NAME, FASTA_ALLELES)
    if list(Path(TEMPDIR, "fasta").glob("insertion*.fasta")):
        # mapping to insertion alleles
        for path_fasta in Path(TEMPDIR, "fasta").glob("insertion*.fasta"):
            name_fasta = path_fasta.stem
            preprocess.mappy_align.output_sam(TEMPDIR, path_fasta, name_fasta, CONTROL, CONTROL_NAME, threads=THREADS)
            preprocess.mappy_align.output_sam(
                TEMPDIR, path_fasta, name_fasta, CONTROL, CONTROL_NAME, preset="splice", threads=THREADS
            )
            preprocess.mappy_align.output_sam(TEMPDIR, path_fasta, name_fasta, SAMPLE, SAMPLE_NAME, threads=THREADS)
            preprocess.mappy_align.output_sam(
                TEMPDIR, path_fasta, name_fasta, SAMPLE, SAMPLE_NAME, preset="splice", threads=THREADS
            )
        # add insertions to FASTA_ALLELES
        for path_fasta in Path(TEMPDIR, "fasta").glob("insertion*.fasta"):
            allele, seq = Path(path_fasta).read_text().strip().split("\n")
            allele = allele.replace(">", "")
            FASTA_ALLELES[allele] = seq
        # MIDSV conversion
        preprocess.call_midsv(TEMPDIR, FASTA_ALLELES, CONTROL_NAME)
        preprocess.call_midsv(TEMPDIR, FASTA_ALLELES, SAMPLE_NAME)
        # Reculculate mutation loci
        preprocess.process_mutation_loci(TEMPDIR, FASTA_ALLELES, CONTROL_NAME)
        preprocess.extract_knockin_loci(TEMPDIR, SAMPLE_NAME)
        preprocess.extract_mutation_loci(TEMPDIR, FASTA_ALLELES, SAMPLE_NAME, CONTROL_NAME)
    ########################################################################
    # Classify alleles
    ########################################################################
    print(f"{_dtnow()}: Classify {arguments['sample']}...", file=sys.stderr)
    classif_sample = classification.classify_alleles(TEMPDIR, FASTA_ALLELES, SAMPLE_NAME)
    with open(Path(TEMPDIR, "classification", f"{SAMPLE_NAME}.pickle"), "wb") as p:
        pickle.dump(classif_sample, p)
    ########################################################################
    # Clustering
    ########################################################################
    print(f"{_dtnow()}: Clustering {arguments['sample']}...", file=sys.stderr)
    clust_sample = clustering.add_labels(classif_sample, TEMPDIR, SAMPLE_NAME, CONTROL_NAME, THREADS)
    clust_sample = clustering.add_readnum(clust_sample)
    clust_sample = clustering.add_percent(clust_sample)
    clust_sample = clustering.update_labels(clust_sample)
    with open(Path(TEMPDIR, "clustering", f"{SAMPLE_NAME}.pickle"), "wb") as p:
        pickle.dump(clust_sample, p)
    ########################################################################
    # Consensus call
    ########################################################################
    print(f"{_dtnow()}: Consensus calling of {arguments['sample']}...", file=sys.stderr)
    # Downsampling to 1000 reads in each LABEL
    MUTATION_LOCI_LABELS = consensus.extract_mutation_loci_by_labels(
        clust_sample, TEMPDIR, FASTA_ALLELES, CONTROL_NAME, SAMPLE_NAME
    )
    clust_subset_sample = consensus.subset_clust(clust_sample, 1000)
    cons_percentage, cons_sequence = consensus.call_consensus(clust_subset_sample, MUTATION_LOCI_LABELS)
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
    report.report_bam.output_bam_sample(TEMPDIR, RESULT_SAMPLE, SAMPLE_NAME, GENOME, GENOME_COODINATES, THREADS)
    for path_bam_igvjs in Path(TEMPDIR, "cache", ".igvjs").glob(f"{CONTROL_NAME}_control.bam*"):
        shutil.copy(path_bam_igvjs, Path(TEMPDIR, "report", ".igvjs", SAMPLE_NAME))
    # VCF
    # working in progress
    ###########################################################
    # Finish call
    ###########################################################
    print(f"{_dtnow()}: \N{teacup without handle} {arguments['sample']} is finished!", file=sys.stderr)
