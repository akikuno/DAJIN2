from __future__ import annotations

import hashlib
import shutil
from pathlib import Path
import json
import midsv

from DAJIN2.core import classification, clustering, consensus, preprocess, report


def parse_args(arguments: dict):
    SAMPLE = arguments["sample"]
    CONTROL = arguments["control"]
    ALLELE = arguments["allele"]
    NAME = arguments["name"]
    THREADS = arguments["threads"]
    try:
        GENOME = arguments["genome"]
    except KeyError:
        GENOME = None
    return SAMPLE, CONTROL, ALLELE, NAME, THREADS, GENOME


def format_inputs(SAMPLE, CONTROL, ALLELE, NAME, GENOME):
    SAMPLE_NAME = preprocess.format_inputs.extract_basename(SAMPLE)
    CONTROL_NAME = preprocess.format_inputs.extract_basename(CONTROL)
    FASTA_ALLELES = preprocess.format_inputs.dictionize_allele(ALLELE)
    TEMPDIR = Path("DAJINResults", ".tempdir", NAME)
    preprocess.format_inputs.make_directories(TEMPDIR, SAMPLE_NAME, CONTROL_NAME)
    if GENOME:
        GENOME_COODINATES = json.loads(Path(TEMPDIR, "cache", "genome_coodinates.jsonl").read_text())
        CHROME_SIZE = int(Path(TEMPDIR, "cache", "chrome_size.txt").read_text())
    else:
        GENOME_COODINATES = None
        CHROME_SIZE = None
    return SAMPLE_NAME, CONTROL_NAME, FASTA_ALLELES, TEMPDIR, GENOME_COODINATES, CHROME_SIZE


def manage_inputs(arguments: dict):
    SAMPLE, CONTROL, ALLELE, NAME, THREADS, GENOME = parse_args(arguments)
    preprocess.validate_inputs.check_files(SAMPLE, CONTROL, ALLELE)
    if GENOME:
        UCSC_URL, GOLDENPATH_URL = preprocess.validate_inputs.check_and_fetch_genome(GENOME)
    # Format inputs -----------------------------
    SAMPLE: str = preprocess.format_inputs.convert_to_posix_path(SAMPLE)
    CONTROL: str = preprocess.format_inputs.convert_to_posix_path(CONTROL)
    ALLELE: str = preprocess.format_inputs.convert_to_posix_path(ALLELE)
    SAMPLE_NAME: str = preprocess.format_inputs.extract_basename(SAMPLE)
    CONTROL_NAME: str = preprocess.format_inputs.extract_basename(CONTROL)
    FASTA_ALLELES: dict = preprocess.format_inputs.dictionize_allele(ALLELE)
    TEMPDIR: Path = Path("DAJINResults", ".tempdir", NAME)
    preprocess.format_inputs.make_directories(TEMPDIR, SAMPLE_NAME, CONTROL_NAME)
    if GENOME:
        GENOME_COODINATES: dict = preprocess.format_inputs.fetch_coodinate(GENOME, UCSC_URL, FASTA_ALLELES["control"])
        CHROME_SIZE: int = preprocess.format_inputs.fetch_chrom_size(GENOME_COODINATES["chr"], GENOME, GOLDENPATH_URL)
        preprocess.format_inputs.cache_coodinates_and_chromsize(TEMPDIR, GENOME, GENOME_COODINATES, CHROME_SIZE)


def execute_control(arguments: dict):
    SAMPLE, CONTROL, ALLELE, NAME, THREADS, GENOME = parse_args(arguments)
    preprocess.validate_inputs.check_files(SAMPLE, CONTROL, ALLELE)
    _, CONTROL_NAME, FASTA_ALLELES, TEMPDIR, GENOME_COODINATES, CHROME_SIZE = format_inputs(
        SAMPLE, CONTROL, ALLELE, NAME, GENOME
    )
    ##########################################################
    # Mapping with minimap2/mappy
    ##########################################################
    for path_fasta in Path(TEMPDIR, "fasta").glob("*.fasta"):
        name_fasta = path_fasta.stem
        preprocess.mappy_align.output_sam(TEMPDIR, path_fasta, name_fasta, CONTROL, CONTROL_NAME)
    if not GENOME:
        path_control = Path(TEMPDIR, "fasta", "control.fasta")
        faidx_control = preprocess.mappy_align.make_faidx(path_control)
        shutil.copy(path_control, Path(TEMPDIR, "report", ".igvjs"))
        Path(TEMPDIR, "report", ".igvjs", "control.fasta.fai").write_text(faidx_control)
    ########################################################################
    # MIDSV conversion
    ########################################################################
    for path_sam in Path(TEMPDIR, "sam").glob(f"{CONTROL_NAME}*"):
        preprocess.call_midsv.output_midsv(TEMPDIR, path_sam, FASTA_ALLELES)
    ########################################################################
    # Classify alleles
    ########################################################################
    path_midsv = Path(TEMPDIR, "midsv").glob(f"{CONTROL_NAME}*")
    classif_control = classification.classify_alleles(path_midsv, CONTROL_NAME)
    ########################################################################
    # Detect Structural variants
    ########################################################################
    for classif in classif_control:
        classif["SV"] = classification.detect_sv(classif["CSSPLIT"], threshold=50)
    ########################################################################
    # Output classification results to cache
    ########################################################################
    midsv.write_jsonl(classif_control, Path(TEMPDIR, "cache", f"{CONTROL_NAME}.jsonl"))
    ########################################################################
    # Output BAM to cache
    ########################################################################
    report.report_bam.output_bam_control(TEMPDIR, CONTROL_NAME, GENOME, GENOME_COODINATES, CHROME_SIZE, THREADS)
    ###############################################################################
    # Cashe inputs (control)
    ###############################################################################
    control_hash = Path(CONTROL).read_bytes()
    control_hash = hashlib.sha256(control_hash).hexdigest()
    PATH_CACHE_HASH = Path(TEMPDIR, "cache", "control_hash.txt")
    PATH_CACHE_HASH.write_text(str(control_hash))


def execute_sample(arguments: dict):
    SAMPLE, CONTROL, ALLELE, NAME, THREADS, GENOME = parse_args(arguments)

    preprocess.validate_inputs.check_files(SAMPLE, CONTROL, ALLELE)

    SAMPLE_NAME, CONTROL_NAME, FASTA_ALLELES, TEMPDIR, GENOME_COODINATES, CHROME_SIZE = format_inputs(
        SAMPLE, CONTROL, ALLELE, NAME, GENOME
    )

    ###############################################################################
    # Mapping with minimap2/mappy
    ###############################################################################

    for path_fasta in Path(TEMPDIR, "fasta").glob("*.fasta"):
        name_fasta = path_fasta.stem
        preprocess.mappy_align.output_sam(TEMPDIR, path_fasta, name_fasta, SAMPLE, SAMPLE_NAME)

    ########################################################################
    # MIDSV conversion
    ########################################################################

    for path_sam in Path(TEMPDIR, "sam").glob(f"{SAMPLE_NAME}*"):
        preprocess.call_midsv.output_midsv(TEMPDIR, path_sam, FASTA_ALLELES)

    ########################################################################
    # Classify alleles
    ########################################################################

    path_midsv = Path(TEMPDIR, "midsv").glob(f"{SAMPLE_NAME}*")
    classif_sample = classification.classify_alleles(path_midsv, SAMPLE_NAME)

    ########################################################################
    # Detect Structural variants
    ########################################################################

    for classif in classif_sample:
        classif["SV"] = classification.detect_sv(classif["CSSPLIT"], threshold=50)

    ########################################################################
    # Clustering
    ########################################################################

    KNOCKIN_LOCI = clustering.find_knockin_loci(TEMPDIR, FASTA_ALLELES, CONTROL_NAME)

    DIFFLOCI_ALLELES, REPETITIVE_DELLOCI = clustering.extract_different_loci(
        TEMPDIR, classif_sample, KNOCKIN_LOCI, FASTA_ALLELES, CONTROL_NAME
    )

    # classif_control = midsv.read_jsonl(Path(TEMPDIR, "cache", f"{CONTROL_NAME}.jsonl"))

    clust_sample = clustering.add_labels(classif_sample, DIFFLOCI_ALLELES)
    clust_sample = clustering.add_readnum(clust_sample)
    clust_sample = clustering.add_percent(clust_sample)
    clust_sample = clustering.update_labels(clust_sample)

    ########################################################################
    # Consensus call
    ########################################################################

    RESULT_SAMPLE, cons_percentage, cons_sequence = consensus.call(
        clust_sample, DIFFLOCI_ALLELES, REPETITIVE_DELLOCI, FASTA_ALLELES
    )

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
