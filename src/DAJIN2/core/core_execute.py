from __future__ import annotations
import hashlib
import shutil
from pathlib import Path
import json
import midsv
from DAJIN2.core import classification, clustering, consensus, preprocess, report


def parse_args(arguments: dict):
    SAMPLE: str = arguments["sample"]
    CONTROL: str = arguments["control"]
    ALLELE: str = arguments["allele"]
    NAME: str = arguments["name"]
    THREADS: int = arguments["threads"]
    try:
        GENOME: str = arguments["genome"]
    except KeyError:
        GENOME = None
    return SAMPLE, CONTROL, ALLELE, NAME, THREADS, GENOME


def format_inputs(SAMPLE, CONTROL, ALLELE, NAME, GENOME, THREADS):
    SAMPLE = preprocess.format_inputs.convert_to_posix_path(SAMPLE)
    CONTROL = preprocess.format_inputs.convert_to_posix_path(CONTROL)
    ALLELE = preprocess.format_inputs.convert_to_posix_path(ALLELE)
    SAMPLE_NAME: str = preprocess.format_inputs.extract_basename(SAMPLE)
    CONTROL_NAME: str = preprocess.format_inputs.extract_basename(CONTROL)
    FASTA_ALLELES: dict = preprocess.format_inputs.dictionize_allele(ALLELE)
    THREADS: int = preprocess.format_inputs.update_threads(THREADS)
    TEMPDIR = Path("DAJINResults", ".tempdir", NAME)
    preprocess.format_inputs.make_directories(TEMPDIR, SAMPLE_NAME, CONTROL_NAME)
    IS_CACHE_CONTROL = preprocess.validate_inputs.exists_cached_control(CONTROL, TEMPDIR)
    IS_CACHE_GENOME = preprocess.validate_inputs.exists_cached_genome(GENOME, TEMPDIR, IS_CACHE_CONTROL)
    UCSC_URL, GOLDENPATH_URL = None, None
    if GENOME:
        if not IS_CACHE_GENOME:
            UCSC_URL, GOLDENPATH_URL = preprocess.validate_inputs.check_and_fetch_genome(GENOME)
            GENOME_COODINATES = preprocess.format_inputs.fetch_coodinate(GENOME, UCSC_URL, FASTA_ALLELES["control"])
            CHROME_SIZE = preprocess.format_inputs.fetch_chrom_size(GENOME_COODINATES["chr"], GENOME, GOLDENPATH_URL)
            preprocess.format_inputs.cache_coodinates_and_chromsize(TEMPDIR, GENOME, GENOME_COODINATES, CHROME_SIZE)
        else:
            GENOME_COODINATES = json.loads(Path(TEMPDIR, "cache", "genome_coodinates.jsonl").read_text())
            CHROME_SIZE = int(Path(TEMPDIR, "cache", "chrome_size.txt").read_text())
    return SAMPLE_NAME, CONTROL_NAME, FASTA_ALLELES, TEMPDIR, GENOME_COODINATES, CHROME_SIZE, THREADS


def execute_control(arguments: dict):
    ###########################################################
    # Preprocess
    ###########################################################
    SAMPLE, CONTROL, ALLELE, NAME, THREADS, GENOME = parse_args(arguments)
    preprocess.validate_inputs.check_files(SAMPLE, CONTROL, ALLELE)
    _, CONTROL_NAME, FASTA_ALLELES, TEMPDIR, GENOME_COODINATES, CHROME_SIZE, THREADS = format_inputs(
        SAMPLE, CONTROL, ALLELE, NAME, GENOME, THREADS
    )
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
    # ============================================================
    # MIDSV conversion
    # ============================================================
    for allele in FASTA_ALLELES:
        preprocess.call_midsv(TEMPDIR, CONTROL_NAME, allele)
    # for path_sam in Path(TEMPDIR, "sam").glob(f"{CONTROL_NAME}_splice_*"):
    #     preprocess.call_midsv(TEMPDIR, path_sam)
    # ============================================================
    # Convert any `N` as deletions other than consecutive `N` from both ends
    # ============================================================
    preprocess.replace_NtoD(TEMPDIR, FASTA_ALLELES, CONTROL_NAME)
    # ============================================================
    # Output MIDSV and BAM to cache
    # ============================================================
    shutil.copytree(Path(TEMPDIR, "midsv"), Path(TEMPDIR, "midsv_control"), dirs_exist_ok=True)
    report.report_bam.output_bam_control(TEMPDIR, CONTROL_NAME, GENOME, GENOME_COODINATES, CHROME_SIZE, THREADS)


def execute_sample(arguments: dict):
    ###########################################################
    # Preprocess
    ###########################################################
    SAMPLE, CONTROL, ALLELE, NAME, THREADS, GENOME = parse_args(arguments)
    preprocess.validate_inputs.check_files(SAMPLE, CONTROL, ALLELE)
    SAMPLE_NAME, CONTROL_NAME, FASTA_ALLELES, TEMPDIR, GENOME_COODINATES, CHROME_SIZE, THREADS = format_inputs(
        SAMPLE, CONTROL, ALLELE, NAME, GENOME, THREADS
    )
    # ============================================================
    # Mapping with minimap2/mappy
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
    shutil.copytree(Path(TEMPDIR, "midsv_control"), Path(TEMPDIR, "midsv"), dirs_exist_ok=True)
    for allele in FASTA_ALLELES:
        preprocess.call_midsv(TEMPDIR, SAMPLE_NAME, allele)
    # for path_sam in Path(TEMPDIR, "sam").glob(f"{SAMPLE_NAME}_splice_*"):
    #     preprocess.call_midsv(TEMPDIR, path_sam)
    # ============================================================
    # Convert any `N` as deletions other than consecutive `N` from both ends
    # ============================================================
    preprocess.replace_NtoD(TEMPDIR, FASTA_ALLELES, SAMPLE_NAME)
    # ============================================================
    # CSSPLITS Error Correction
    # ============================================================
    MUTATION_LOCI_ALLELES = preprocess.extract_mutation_loci(TEMPDIR, FASTA_ALLELES, CONTROL_NAME, SAMPLE_NAME)
    preprocess.correct_sequence_error(TEMPDIR, FASTA_ALLELES, CONTROL_NAME, SAMPLE_NAME, MUTATION_LOCI_ALLELES)
    KNOCKIN_LOCI_ALLELES = preprocess.extract_knockin_loci(TEMPDIR)
    # preprocess.correct_knockin.execute(TEMPDIR, FASTA_ALLELES, CONTROL_NAME, SAMPLE_NAME)
    ########################################################################
    # Classify alleles
    ########################################################################
    print("Classify...")
    classif_sample = classification.classify_alleles(TEMPDIR, SAMPLE_NAME)
    # for classif in classif_sample:
    #     classif["SV"] = classification.detect_sv(classif["CSSPLIT"], threshold=50)
    ########################################################################
    # Clustering
    ########################################################################
    print("Clustering...")
    clust_sample = clustering.add_labels(classif_sample, TEMPDIR, CONTROL_NAME, MUTATION_LOCI_ALLELES, KNOCKIN_LOCI_ALLELES, THREADS)
    clust_sample = clustering.add_readnum(clust_sample)
    clust_sample = clustering.add_percent(clust_sample)
    clust_sample = clustering.update_labels(clust_sample)
    ########################################################################
    # Consensus call
    ########################################################################
    print("Consensus call...")
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
