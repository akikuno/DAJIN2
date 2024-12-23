from __future__ import annotations

import logging
import shutil
from pathlib import Path

from DAJIN2.core import classification, clustering, consensus, preprocess, report
from DAJIN2.core.preprocess.input_formatter import FormattedInputs
from DAJIN2.utils import fastx_handler, io

logger = logging.getLogger(__name__)


###########################################################
# main
###########################################################


def execute_control(arguments: dict):
    logger.info(f"{arguments['control']} is now processing...")

    ###########################################################
    # Setup arguments
    ###########################################################
    ARGS: FormattedInputs = preprocess.format_inputs(arguments)
    preprocess.create_temporal_directories(ARGS.tempdir, ARGS.control_name, is_control=True)
    preprocess.create_report_directories(ARGS.tempdir, ARGS.control_name, is_control=True)
    io.cache_control_hash(ARGS.tempdir, ARGS.path_allele)

    ###########################################################
    # Check caches
    ###########################################################
    if Path(ARGS.tempdir, "report", "BAM", ARGS.control_name, f"{ARGS.control_name}.bam").exists():
        logger.info(f"{arguments['control']} is already preprocessed and reuse the results for the current run...")
        return None

    logger.info(f"Preprocess {arguments['control']}...")

    ###########################################################
    # Concatenate fastq files
    ###########################################################

    fastx_handler.save_inputs_as_single_fastq(ARGS, is_control=True)

    path_fastq = Path(ARGS.tempdir, ARGS.control_name, "fastq", f"{ARGS.control_name}.fastq.gz")
    fastx_handler.overwrite_with_downsampled_fastq(path_fastq, num_reads=20_000)

    ###########################################################
    # Export fasta files as single-FASTA format
    ###########################################################
    fastx_handler.export_fasta_files(ARGS, is_control=True)

    ###########################################################
    # Separate fastq files by sequence error
    ###########################################################
    paths_fasta = Path(ARGS.tempdir, ARGS.control_name, "fasta").glob("control.fasta")
    preprocess.generate_sam(ARGS, paths_fasta, is_control=True, is_sv=False)
    preprocess.generate_midsv(ARGS, is_control=True, is_sv=False)
    # ============================================================
    # Detect sequence error reads
    # ============================================================
    preprocess.detect_sequence_error_reads(ARGS, is_control=True)
    preprocess.split_fastq_by_sequence_error(ARGS, is_control=True)

    # ============================================================
    # Save subsetted fastq if the read number is too large (> 10,000 reads)
    # ============================================================
    path_fastq = Path(ARGS.tempdir, ARGS.control_name, "fastq", f"{ARGS.control_name}.fastq.gz")
    path_fastq_with_error = Path(
        ARGS.tempdir, ARGS.control_name, "fastq", f"{ARGS.control_name}_sequence_error.fastq.gz"
    )
    fastx_handler.overwrite_with_downsampled_fastq(path_fastq, num_reads=10_000)
    fastx_handler.overwrite_with_downsampled_fastq(path_fastq_with_error, num_reads=10_000)

    ###########################################################
    # Mapping filtered control reads
    ###########################################################

    # ============================================================
    # Mapping using mappy
    # ============================================================
    paths_fasta = Path(ARGS.tempdir, ARGS.control_name, "fasta").glob("*.fasta")
    preprocess.generate_sam(ARGS, paths_fasta, is_control=True, is_sv=False)

    # ============================================================
    # MIDSV conversion
    # ============================================================
    preprocess.generate_midsv(ARGS, is_control=True, is_sv=False)

    ###########################################################
    # Prepare data to `extract mutaion loci`
    ###########################################################
    preprocess.cache_mutation_loci(ARGS, is_control=True)

    ###########################################################
    # Output BAM files
    ###########################################################
    logger.info(f"Output BAM files of {arguments['control']}...")
    report.bam_exporter.export_to_bam(
        ARGS.tempdir, ARGS.control_name, ARGS.genome_coordinates, ARGS.threads, is_control=True
    )
    ###########################################################
    # Finish call
    ###########################################################
    logger.info(f"\N{TEACUP WITHOUT HANDLE} {arguments['control']} is finished!")


def execute_sample(arguments: dict):
    logger.info(f"{arguments['sample']} is now processing...")

    ###########################################################
    # Preprocess
    ###########################################################

    ARGS: FormattedInputs = preprocess.format_inputs(arguments)
    preprocess.create_temporal_directories(ARGS.tempdir, ARGS.sample_name, is_control=False)
    preprocess.create_report_directories(ARGS.tempdir, ARGS.sample_name, is_control=False)

    logger.info(f"Preprocess {arguments['sample']}...")

    # ============================================================
    # Merge fastq files
    # ============================================================

    fastx_handler.save_inputs_as_single_fastq(ARGS, is_control=False)

    # Save subsetted fastq if the read number is too large (> 10,000 reads)
    path_fastq = Path(ARGS.tempdir, ARGS.sample_name, "fastq", f"{ARGS.sample_name}.fastq.gz")
    fastx_handler.overwrite_with_downsampled_fastq(path_fastq, num_reads=10_000)

    # ============================================================
    # Mapping with mappy
    # ============================================================

    for path_fasta in Path(ARGS.tempdir, ARGS.control_name, "fasta").glob("*.fasta"):
        shutil.copy(path_fasta, Path(ARGS.tempdir, ARGS.sample_name, "fasta"))

    paths_fasta = Path(ARGS.tempdir, ARGS.sample_name, "fasta").glob("*.fasta")
    preprocess.generate_sam(ARGS, paths_fasta, is_control=False, is_sv=False)

    # ============================================================
    # MIDSV conversion
    # ============================================================

    preprocess.generate_midsv(ARGS, is_control=False, is_sv=False)

    ###########################################################
    # Remove sequence error reads from MIDSV files
    ###########################################################
    preprocess.detect_sequence_error_reads(ARGS, is_control=False)
    preprocess.split_fastq_by_sequence_error(ARGS, is_control=False)
    preprocess.replace_midsv_without_sequence_errors(ARGS)

    # ============================================================
    # Extract mutation loci
    # ============================================================
    preprocess.extract_knockin_loci(ARGS.tempdir, ARGS.sample_name)
    preprocess.cache_mutation_loci(ARGS, is_control=False)

    # ============================================================
    # Detect and mapping SV alleles
    # ============================================================
    paths_predefined_fasta: set[str] = {
        str(Path(ARGS.tempdir, ARGS.sample_name, "fasta", f"{allele}.fasta")) for allele in ARGS.fasta_alleles.keys()
    }

    preprocess.detect_insertions(ARGS.tempdir, ARGS.sample_name, ARGS.control_name, ARGS.fasta_alleles)
    # preprocess.detect_inversions(ARGS.tempdir, ARGS.sample_name, ARGS.control_name, ARGS.fasta_alleles) #TODO

    paths_sv_fasta = {str(p) for p in Path(ARGS.tempdir, ARGS.sample_name, "fasta").glob("insertion*.fasta")}
    # paths_sv_fasta |= {str(p) for p in Path(ARGS.tempdir, ARGS.sample_name, "fasta").glob("inversion*.fasta")}

    paths_sv_fasta -= paths_predefined_fasta

    if paths_sv_fasta:
        # mapping to SV alleles
        preprocess.generate_sam(ARGS, paths_sv_fasta, is_control=True, is_sv=True)
        preprocess.generate_sam(ARGS, paths_sv_fasta, is_control=False, is_sv=True)
        # add SV alleles to ARGS.fasta_alleles
        for path_fasta_sv in paths_sv_fasta:
            allele, seq = Path(path_fasta_sv).read_text().strip().split("\n")
            allele = allele.replace(">", "")
            ARGS.fasta_alleles[allele] = seq
        # MIDSV conversion
        preprocess.generate_midsv(ARGS, is_control=True, is_sv=True)
        preprocess.generate_midsv(ARGS, is_control=False, is_sv=True)
        # Reculculate mutation loci
        preprocess.cache_mutation_loci(ARGS, is_control=True)
        preprocess.extract_knockin_loci(ARGS.tempdir, ARGS.sample_name)
        preprocess.cache_mutation_loci(ARGS, is_control=False)

    io.save_pickle(ARGS.fasta_alleles, Path(ARGS.tempdir, ARGS.sample_name, "fasta", "fasta_alleles.pickle"))

    ########################################################################
    # Classify alleles
    ########################################################################

    logger.info(f"Classify {arguments['sample']}...")

    classif_sample = classification.classify_alleles(ARGS.tempdir, ARGS.fasta_alleles, ARGS.sample_name)

    io.save_pickle(classif_sample, Path(ARGS.tempdir, ARGS.sample_name, "classification", "classif_sample.pickle"))

    ########################################################################
    # Clustering
    ########################################################################

    logger.info(f"Clustering {arguments['sample']}...")

    labels = clustering.extract_labels(classif_sample, ARGS.tempdir, ARGS.sample_name, ARGS.control_name)
    clust_sample = clustering.add_labels(classif_sample, labels)
    clust_sample = clustering.add_readnum(clust_sample)
    clust_sample = clustering.add_percent(clust_sample)
    clust_sample = clustering.update_labels(clust_sample)

    io.save_pickle(clust_sample, Path(ARGS.tempdir, ARGS.sample_name, "clustering", "clust_sample.pickle"))

    ########################################################################
    # Consensus call
    ########################################################################

    logger.info(f"Consensus calling of {arguments['sample']}...")

    # Remove minor alleles with fewer than 5 reads or less than 0.5%
    clust_sample_removed = consensus.remove_minor_alleles(clust_sample)

    # Downsampling to 1000 reads in each LABEL
    clust_downsampled = consensus.downsample_by_label(clust_sample_removed, 1000)

    consensus.cache_mutation_loci(ARGS, clust_downsampled)

    cons_percentage, cons_sequence = consensus.call_consensus(
        ARGS.tempdir, ARGS.sample_name, ARGS.fasta_alleles, clust_downsampled
    )

    allele_names = consensus.call_allele_name(cons_sequence, cons_percentage, ARGS.fasta_alleles)
    cons_percentage = consensus.update_key_by_allele_name(cons_percentage, allele_names)
    cons_sequence = consensus.update_key_by_allele_name(cons_sequence, allele_names)

    RESULT_SAMPLE = consensus.add_key_by_allele_name(clust_sample_removed, allele_names)
    RESULT_SAMPLE.sort(key=lambda x: x["LABEL"])

    io.save_pickle(cons_percentage, Path(ARGS.tempdir, ARGS.sample_name, "consensus", "cons_percentage.pickle"))
    io.save_pickle(cons_sequence, Path(ARGS.tempdir, ARGS.sample_name, "consensus", "cons_sequence.pickle"))

    ########################################################################
    # Output Report：RESULT/FASTA/HTML/BAM
    ########################################################################

    logger.info(f"Output reports of {arguments['sample']}...")

    # RESULT
    io.write_jsonl(RESULT_SAMPLE, Path(ARGS.tempdir, "result", f"{ARGS.sample_name}.jsonl"))
    # FASTA
    report.sequence_exporter.export_to_fasta(ARGS.tempdir, ARGS.sample_name, cons_sequence)
    report.sequence_exporter.export_reference_to_fasta(ARGS.tempdir, ARGS.sample_name)
    # HTML
    report.sequence_exporter.export_to_html(ARGS.tempdir, ARGS.sample_name, cons_percentage)
    # CSV (Allele Info)
    report.mutation_exporter.export_to_csv(ARGS.tempdir, ARGS.sample_name, ARGS.genome_coordinates, cons_percentage)
    # BAM
    report.bam_exporter.export_to_bam(
        ARGS.tempdir, ARGS.sample_name, ARGS.genome_coordinates, ARGS.threads, RESULT_SAMPLE
    )
    for path_bam_igvjs in Path(ARGS.tempdir, "cache", ".igvjs").glob(f"{ARGS.control_name}_control.bam*"):
        shutil.copy(path_bam_igvjs, Path(ARGS.tempdir, "report", ".igvjs", ARGS.sample_name))
    # VCF
    # working in progress

    ###########################################################
    # Finish call
    ###########################################################

    logger.info(f"\N{TEACUP WITHOUT HANDLE} {arguments['sample']} is finished!")
