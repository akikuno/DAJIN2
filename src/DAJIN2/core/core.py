from __future__ import annotations
import hashlib
import shutil
import midsv
from pathlib import Path
from . import preprocess
from . import classification
from . import clustering
from . import consensus
from . import report


def execute(arguments: dict) -> None:
    SAMPLE = arguments["sample"]
    CONTROL = arguments["control"]
    ALLELE = arguments["allele"]
    NAME = arguments["name"]
    THREADS = arguments["threads"]
    try:
        GENOME = arguments["genome"]
    except KeyError:
        GENOME = False
        GENOME_COODINATES = None
        CHROME_SIZE = None

    ##########################################################
    # Check inputs
    ##########################################################
    preprocess.check_inputs.check_files(SAMPLE, CONTROL, ALLELE)
    TEMPDIR = Path("DAJINResults", ".tempdir", NAME)
    EXISTS_CACHED_CONTROL = preprocess.check_inputs.exists_cached_control(CONTROL, TEMPDIR)
    EXISTS_CACHED_GENOME = preprocess.check_inputs.exists_cached_genome(GENOME, TEMPDIR, EXISTS_CACHED_CONTROL)
    UCSC_URL, GOLDENPATH_URL = None, None
    if GENOME and not EXISTS_CACHED_GENOME:
        UCSC_URL, GOLDENPATH_URL = preprocess.check_inputs.check_and_fetch_genome(GENOME)

    ##########################################################
    # Format inputs
    ##########################################################
    SAMPLE_NAME = preprocess.format_inputs.extract_basename(SAMPLE)
    CONTROL_NAME = preprocess.format_inputs.extract_basename(CONTROL)
    DICT_ALLELE = preprocess.format_inputs.dictionize_allele(ALLELE)

    preprocess.format_inputs.make_directories(TEMPDIR)

    if GENOME:
        GENOME_COODINATES, CHROME_SIZE = preprocess.format_inputs.get_coodinates_and_chromsize(
            TEMPDIR, GENOME, DICT_ALLELE, UCSC_URL, GOLDENPATH_URL, EXISTS_CACHED_GENOME
        )
        if not EXISTS_CACHED_GENOME:
            preprocess.format_inputs.cache_coodinates_and_chromsize(TEMPDIR, GENOME, GENOME_COODINATES, CHROME_SIZE)

    ################################################################################
    # Export fasta files as single-FASTA format
    ################################################################################
    # TODO: use yeild, not export
    for identifier, sequence in DICT_ALLELE.items():
        contents = "\n".join([">" + identifier, sequence]) + "\n"
        output_fasta = Path(TEMPDIR, "fasta", f"{identifier}.fasta")
        output_fasta.write_text(contents)

    ###############################################################################
    # Mapping with minimap2/mappy
    ###############################################################################
    for path_fasta in Path(TEMPDIR, "fasta").glob("*.fasta"):
        name_fasta = path_fasta.stem
        if name_fasta not in set(DICT_ALLELE.keys()):
            continue
        if not EXISTS_CACHED_CONTROL:
            preprocess.mappy_align.output_sam(TEMPDIR, path_fasta, name_fasta, CONTROL, CONTROL_NAME)
        preprocess.mappy_align.output_sam(TEMPDIR, path_fasta, name_fasta, SAMPLE, SAMPLE_NAME)

    if not GENOME:
        path_control = Path(TEMPDIR, "fasta", "control.fasta")
        faidx_control = preprocess.mappy_align.make_faidx(path_control)
        shutil.copy(path_control, Path(TEMPDIR, "report", ".igvjs"))
        Path(TEMPDIR, "report", ".igvjs", "control.fasta.fai").write_text(faidx_control)

    ########################################################################
    # MIDSV conversion
    ########################################################################
    if not EXISTS_CACHED_CONTROL:
        for path_sam in Path(TEMPDIR, "sam").glob(f"{CONTROL_NAME}*"):
            preprocess.calc_midsv.output_midsv(TEMPDIR, path_sam, DICT_ALLELE)

    for path_sam in Path(TEMPDIR, "sam").glob(f"{SAMPLE_NAME}*"):
        preprocess.calc_midsv.output_midsv(TEMPDIR, path_sam, DICT_ALLELE)

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

    MASKS_CONTROL = clustering.mask_control(TEMPDIR, DICT_ALLELE, CONTROL_NAME)

    DIFFLOCI_ALLELES, REPETITIVE_DELLOCI = clustering.extract_different_loci(
        TEMPDIR, classif_sample, MASKS_CONTROL, DICT_ALLELE, CONTROL_NAME
    )
    clust_sample = clustering.add_labels(classif_sample, DIFFLOCI_ALLELES)
    clust_sample = clustering.add_readnum(clust_sample)
    clust_sample = clustering.add_percent(clust_sample)
    clust_sample = clustering.update_labels(clust_sample)

    ########################################################################
    # Consensus call
    ########################################################################

    RESULT_SAMPLE, cons_percentage, cons_sequence = consensus.call(
        clust_sample, DIFFLOCI_ALLELES, REPETITIVE_DELLOCI, DICT_ALLELE
    )
    ########################################################################
    # Output Report：RESULT/FASTA/HTML/BAM/VCF
    ########################################################################
    # RESULT
    midsv.write_jsonl(RESULT_SAMPLE, Path(TEMPDIR, "result", f"{SAMPLE_NAME}.jsonl"))

    # FASTA
    for header, cons_seq in cons_sequence.items():
        cons_fasta = report.report_files.to_fasta(header, cons_seq)
        Path(TEMPDIR, "report", "FASTA", f"{SAMPLE_NAME}_{header}.fasta").write_text(cons_fasta)

    # HTML
    for header, cons_per in cons_percentage.items():
        cons_html = report.report_files.to_html(header, cons_per)
        Path(TEMPDIR, "report", "HTML", f"{SAMPLE_NAME}_{header}.html").write_text(cons_html)

    # BAM
    if not EXISTS_CACHED_CONTROL:
        report.report_bam.output_bam_control(TEMPDIR, CONTROL_NAME, GENOME, GENOME_COODINATES, CHROME_SIZE, THREADS)

    report.report_bam.output_bam_sample(
        TEMPDIR, RESULT_SAMPLE, SAMPLE_NAME, GENOME, GENOME_COODINATES, CHROME_SIZE, THREADS
    )

    # VCF
    # working in progress

    ###############################################################################
    # Cashe inputs (control)
    ###############################################################################

    if not EXISTS_CACHED_CONTROL:
        control_hash = Path(CONTROL).read_bytes()
        control_hash = hashlib.sha256(control_hash).hexdigest()
        PATH_CACHE_HASH = Path(TEMPDIR, "cache", "control_hash.txt")
        PATH_CACHE_HASH.write_text(str(control_hash))

