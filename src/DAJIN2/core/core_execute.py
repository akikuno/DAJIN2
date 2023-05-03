from __future__ import annotations

import json
# limit max memory usage
import os
import pickle
import resource
import shutil
import sys
from datetime import datetime
from pathlib import Path

import midsv

from DAJIN2.core import (classification, clustering, consensus, preprocess,
                         report)

mem_bytes = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')
resource.setrlimit(resource.RLIMIT_DATA, (int(mem_bytes * 9/10), -1))


def _parse_arguments(arguments: dict):
    SAMPLE: str = arguments["sample"]
    CONTROL: str = arguments["control"]
    ALLELE: str = arguments["allele"]
    NAME: str = arguments["name"]
    THREADS: int = arguments["threads"]
    if "genome" in arguments:
        GENOME: str = arguments["genome"]
        URL_UCSC = arguments["ucsc"]
        URL_GOLDENPATH = arguments["goldenpath"]
    else:
        GENOME = ""
        URL_UCSC = ""
        URL_GOLDENPATH = ""
    return SAMPLE, CONTROL, ALLELE, NAME, THREADS, GENOME, URL_UCSC, URL_GOLDENPATH


def _format_inputs(arguments: dict):
    SAMPLE, CONTROL, ALLELE, NAME, THREADS, GENOME, URL_UCSC, URL_GOLDENPATH = _parse_arguments(arguments)
    SAMPLE = preprocess.format_inputs.convert_to_posix_path(SAMPLE)
    CONTROL = preprocess.format_inputs.convert_to_posix_path(CONTROL)
    ALLELE = preprocess.format_inputs.convert_to_posix_path(ALLELE)

    SAMPLE_NAME: str = preprocess.format_inputs.extract_basename(SAMPLE)
    CONTROL_NAME: str = preprocess.format_inputs.extract_basename(CONTROL)
    FASTA_ALLELES: dict = preprocess.format_inputs.dictionize_allele(ALLELE)

    TEMPDIR = Path("DAJINResults", ".tempdir", NAME)
    SUBDIRS = ["cache", "fasta", "sam", "midsv", "report", "result", "mutation_loci"]
    preprocess.format_inputs.make_directories(TEMPDIR, SUBDIRS, SAMPLE_NAME, CONTROL_NAME)

    IS_CACHE_CONTROL = preprocess.validate_inputs.exists_cached_control(CONTROL, TEMPDIR)
    IS_CACHE_GENOME = preprocess.validate_inputs.exists_cached_genome(GENOME, TEMPDIR, IS_CACHE_CONTROL)
    if GENOME:
        if not IS_CACHE_GENOME:
            GENOME_COODINATES = preprocess.format_inputs.fetch_coodinate(GENOME, URL_UCSC, FASTA_ALLELES["control"])
            CHROME_SIZE = preprocess.format_inputs.fetch_chrom_size(GENOME_COODINATES["chr"], GENOME, URL_GOLDENPATH)
            preprocess.format_inputs.cache_coodinates_and_chromsize(TEMPDIR, GENOME, GENOME_COODINATES, CHROME_SIZE)
        else:
            GENOME_COODINATES = json.loads(Path(TEMPDIR, "cache", "genome_coodinates.jsonl").read_text())
            CHROME_SIZE = int(Path(TEMPDIR, "cache", "chrome_size.txt").read_text())
    return SAMPLE_NAME, CONTROL_NAME, FASTA_ALLELES, TEMPDIR, GENOME_COODINATES, CHROME_SIZE, THREADS


def _dtnow() -> str:
    return datetime.now().strftime('%Y-%m-%d %H:%M:%S')

def execute_control(arguments: dict):
    print(f"{arguments['control']} is now processing...")
    ###########################################################
    # Preprocess
    ###########################################################
    SAMPLE, CONTROL, ALLELE, NAME, THREADS, GENOME, _, _ = _parse_arguments(arguments)
    # preprocess.validate_inputs.check_files(SAMPLE, CONTROL, ALLELE)
    _, CONTROL_NAME, FASTA_ALLELES, TEMPDIR, GENOME_COODINATES, CHROME_SIZE, THREADS = _format_inputs(
        arguments
    )
    ###########################################################
    # Save Caches
    ###########################################################
    done_midsv = Path(TEMPDIR, "midsv", f"{CONTROL_NAME}.plk").exists()
    done_bam = Path(TEMPDIR, "report", "BAM", CONTROL_NAME, f"{CONTROL_NAME}.bam").exists()
    if done_midsv and done_bam:
        print(f"{arguments['control']} is finished...")
        return
    print(f"{_dtnow()}: Preprocess {CONTROL_NAME}...")
    # ============================================================
    # Export fasta files as single-FASTA format
    # ============================================================
    for identifier, sequence in FASTA_ALLELES.items():
        contents = "\n".join([">" + identifier, sequence]) + "\n"
        output_fasta = Path(TEMPDIR, "fasta", f"{identifier}.fasta")
        output_fasta.write_text(contents)
    print(f"{_dtnow()}: Mapping {CONTROL_NAME}...")
    # ============================================================
    # Mapping using mappy
    # ============================================================
    for path_fasta in Path(TEMPDIR, "fasta").glob("*.fasta"):
        name_fasta = path_fasta.stem
        preprocess.mappy_align.output_sam(TEMPDIR, path_fasta, name_fasta, CONTROL, CONTROL_NAME, threads=THREADS)
        preprocess.mappy_align.output_sam(
            TEMPDIR, path_fasta, name_fasta, CONTROL, CONTROL_NAME, preset="splice", threads=THREADS
        )
    # ============================================================
    # MIDSV conversion
    # ============================================================
    print(f"{_dtnow()}: Call MIDSV {CONTROL_NAME}...")
    preprocess.call_midsv(TEMPDIR, FASTA_ALLELES, CONTROL_NAME)
    ###########################################################
    # Save MIDSV and BAM
    ###########################################################
    # with open(Path(TEMPDIR, "midsv", f"{CONTROL_NAME}.plk"), 'wb') as p:
    #     pickle.dump(midsv_control_alleles, p)
    report.report_bam.output_bam_control(TEMPDIR, CONTROL_NAME, GENOME, GENOME_COODINATES, CHROME_SIZE, THREADS)
    print(f"{arguments['control']} is finished!")


def execute_sample(arguments: dict):
    print(f"{arguments['sample']} is now processing...")
    ###########################################################
    # Preprocess
    ###########################################################
    SAMPLE, CONTROL, ALLELE, NAME, THREADS, GENOME, _, _ = _parse_arguments(arguments)
    # preprocess.validate_inputs.check_files(SAMPLE, CONTROL, ALLELE)
    SAMPLE_NAME, CONTROL_NAME, FASTA_ALLELES, TEMPDIR, GENOME_COODINATES, CHROME_SIZE, THREADS = _format_inputs(
        arguments
    )
    print(f"{_dtnow()}: Preprocess {SAMPLE_NAME}...")
    # ============================================================
    # Mapping with mappy
    # ============================================================
    print(f"{_dtnow()}: Mapping {SAMPLE_NAME}...") #!====================
    for path_fasta in Path(TEMPDIR, "fasta").glob("*.fasta"):
        name_fasta = path_fasta.stem
        preprocess.mappy_align.output_sam(TEMPDIR, path_fasta, name_fasta, SAMPLE, SAMPLE_NAME, threads=THREADS)
        preprocess.mappy_align.output_sam(
            TEMPDIR, path_fasta, name_fasta, SAMPLE, SAMPLE_NAME, preset="splice", threads=THREADS
        )
    # ============================================================
    # MIDSV conversion
    # ============================================================
    print(f"{_dtnow()}: Call MIDSV {SAMPLE_NAME}...") #!====================
    preprocess.call_midsv(TEMPDIR, FASTA_ALLELES, SAMPLE_NAME)
    # with open(Path(TEMPDIR, "midsv", f"{SAMPLE_NAME}.plk"), 'wb') as p:
    #     pickle.dump(midsv_sample_alleles, p)
    # with open(Path(TEMPDIR, "midsv", f"{CONTROL_NAME}.plk"), 'rb') as p:
    #     midsv_control_alleles = pickle.load(p)
    # ============================================================
    # CSSPLITS Error Correction
    # ============================================================
    print(f"{_dtnow()}: extract_mutation_loci {SAMPLE_NAME}...") #!====================
    try:
        MUTATION_LOCI_ALLELES = preprocess.extract_mutation_loci(TEMPDIR, FASTA_ALLELES, SAMPLE_NAME, CONTROL_NAME)
        # MUTATION_LOCI_ALLELES = preprocess.extract_mutation_loci(midsv_sample_alleles, midsv_control_alleles)
    except Exception as e:
        print(f"{SAMPLE_NAME} is failed at extract mutation loci.")
        print(type(e), e)
        exit(1)
    with open(Path(TEMPDIR, "mutation_loci", f"{SAMPLE_NAME}.plk"), 'wb') as p:
        pickle.dump(MUTATION_LOCI_ALLELES, p)
    exit(0) #!====================

    print(f"{_dtnow()}: correct_sequence_error {SAMPLE_NAME}...") #!====================
    try:
        midsv_alleles_corrected = preprocess.correct_sequence_error(midsv_sample_alleles, midsv_control_alleles, FASTA_ALLELES, MUTATION_LOCI_ALLELES)
    except Exception as e:
        print(f"{SAMPLE_NAME} is failed at `correct_sequence_error`.")
        print(type(e), e)
        exit(1)
    print(f"{SAMPLE_NAME} is finished by correct_sequence_error.")
    KNOCKIN_LOCI_ALLELES = preprocess.extract_knockin_loci(TEMPDIR)
    # preprocess.correct_knockin.execute(TEMPDIR, FASTA_ALLELES, CONTROL_NAME, SAMPLE_NAME)
    ########################################################################
    # Classify alleles
    ########################################################################
    print(f"{_dtnow()}: Classify {SAMPLE_NAME}...")
    classif_sample = classification.classify_alleles(midsv_alleles_corrected["sample"], TEMPDIR)
    ########################################################################
    # Clustering
    ########################################################################
    print(f"{_dtnow()}: Clustering {SAMPLE_NAME}...")
    clust_sample = clustering.add_labels(classif_sample, midsv_alleles_corrected["control"], MUTATION_LOCI_ALLELES, KNOCKIN_LOCI_ALLELES, THREADS)
    clust_sample = clustering.add_readnum(clust_sample)
    clust_sample = clustering.add_percent(clust_sample)
    clust_sample = clustering.update_labels(clust_sample)
    ########################################################################
    # Consensus call
    ########################################################################
    print(f"{_dtnow()}: Consensus calling {SAMPLE_NAME}......")
    cons_percentage, cons_sequence = consensus.call_consensus(clust_sample)
    allele_names = consensus.call_allele_name(cons_sequence, cons_percentage, FASTA_ALLELES)
    cons_percentage = consensus.update_key_by_allele_name(cons_percentage, allele_names)
    cons_sequence = consensus.update_key_by_allele_name(cons_sequence, allele_names)
    RESULT_SAMPLE = consensus.add_key_by_allele_name(clust_sample, allele_names)
    RESULT_SAMPLE.sort(key=lambda x: x["LABEL"])
    ########################################################################
    # Output Report：RESULT/FASTA/HTML/BAM/VCF
    ########################################################################
    # RESULT
    midsv.write_jsonl(RESULT_SAMPLE, Path(TEMPDIR, "result", f"{SAMPLE_NAME}.jsonl"))
    # FASTA
    for header, cons_seq in cons_sequence.items():
        cons_fasta = report.report_files.to_fasta(header, cons_seq)
        Path(TEMPDIR, "report", "FASTA", SAMPLE_NAME, f"{SAMPLE_NAME}_{header}.fasta").write_text(cons_fasta)
    # HTML
    for header, cons_per in cons_percentage.items():
        cons_html = report.report_files.to_html(header, cons_per)
        Path(TEMPDIR, "report", "HTML", SAMPLE_NAME, f"{SAMPLE_NAME}_{header}.html").write_text(cons_html)
    # BAM
    report.report_bam.output_bam_sample(
        TEMPDIR, RESULT_SAMPLE, SAMPLE_NAME, GENOME, GENOME_COODINATES, CHROME_SIZE, THREADS
    )
    for path_bam_igvjs in Path(TEMPDIR, "cache", ".igvjs").glob(f"{CONTROL_NAME}_control.bam*"):
        shutil.copy(path_bam_igvjs, Path(TEMPDIR, "report", ".igvjs", SAMPLE_NAME))
    # VCF
    # working in progress
    print(f"{arguments['sample']} is finished...")
