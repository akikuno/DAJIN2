from __future__ import annotations

import logging
import shutil
from pathlib import Path

from DAJIN2.core import classification, clustering, consensus, preprocess, report
from DAJIN2.core.preprocess.infrastructure.input_formatter import FormattedInputs
from DAJIN2.utils import allele_handler, fastx_handler, fileio

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
    allele_handler.save_allele_name_map(ARGS.tempdir, ARGS.fasta_alleles.keys())
    fileio.cache_file_hash(ARGS.path_allele, Path(ARGS.tempdir, "cache", "hash_allele.txt"))

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
    control_allele_key = allele_handler.to_allele_key("control")
    paths_fasta = [Path(ARGS.tempdir, ARGS.control_name, "fasta", f"{control_allele_key}.fasta")]
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
    if ARGS.control_coordinate_scoring:
        control_allele_key = allele_handler.to_allele_key("control")
        paths_fasta = [Path(ARGS.tempdir, ARGS.control_name, "fasta", f"{control_allele_key}.fasta")]
    else:
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
    allele_handler.save_allele_name_map(ARGS.tempdir, ARGS.fasta_alleles.keys())

    logger.info(f"Preprocess {arguments['sample']}...")

    # ============================================================
    # Merge fastq files
    # ============================================================

    fastx_handler.save_inputs_as_single_fastq(ARGS, is_control=False)

    # Save subsetted fastq if the read number is too large (> 100,000 reads)
    path_fastq = Path(ARGS.tempdir, ARGS.sample_name, "fastq", f"{ARGS.sample_name}.fastq.gz")
    fastx_handler.overwrite_with_downsampled_fastq(path_fastq, num_reads=100_000)

    if ARGS.control_coordinate_scoring:
        _execute_sample_with_control_coordinate_scoring(arguments, ARGS)
        return None

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
        str(
            Path(
                ARGS.tempdir,
                ARGS.sample_name,
                "fasta",
                f"{allele_handler.to_allele_key(allele)}.fasta",
            )
        )
        for allele in ARGS.fasta_alleles.keys()
    }

    preprocess.detect_sv_alleles(
        ARGS.tempdir, ARGS.sample_name, ARGS.control_name, ARGS.fasta_alleles, sv_type="insertion"
    )
    preprocess.detect_sv_alleles(
        ARGS.tempdir, ARGS.sample_name, ARGS.control_name, ARGS.fasta_alleles, sv_type="deletion"
    )
    preprocess.detect_sv_alleles(
        ARGS.tempdir, ARGS.sample_name, ARGS.control_name, ARGS.fasta_alleles, sv_type="inversion"
    )

    paths_all_fasta = {str(p) for p in Path(ARGS.tempdir, ARGS.sample_name, "fasta").glob("*.fasta")}
    paths_sv_fasta = paths_all_fasta - paths_predefined_fasta

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

    allele_handler.save_allele_name_map(ARGS.tempdir, ARGS.fasta_alleles.keys())
    fileio.save_pickle(ARGS.fasta_alleles, Path(ARGS.tempdir, ARGS.sample_name, "fasta", "fasta_alleles.pickle"))

    ########################################################################
    # Classify alleles
    ########################################################################

    logger.info(f"Classify {arguments['sample']}...")

    classif_sample = classification.classify_alleles(
        ARGS.tempdir, ARGS.fasta_alleles, ARGS.sample_name, ARGS.no_filter
    )

    _complete_sample_analysis(arguments, ARGS, classif_sample)


def _execute_sample_with_control_coordinate_scoring(arguments: dict, ARGS: FormattedInputs) -> None:
    from DAJIN2.core.classification import control_coordinate

    control_key = allele_handler.to_allele_key("control")
    original_alleles = set(ARGS.fasta_alleles)

    path_control_fasta = Path(ARGS.tempdir, ARGS.control_name, "fasta", f"{control_key}.fasta")
    shutil.copy(path_control_fasta, Path(ARGS.tempdir, ARGS.sample_name, "fasta"))

    ###########################################################
    # Build sample/control MIDSV on the control coordinate.
    ###########################################################
    preprocess.generate_sam(ARGS, [Path(ARGS.tempdir, ARGS.sample_name, "fasta", f"{control_key}.fasta")])
    preprocess.generate_midsv(ARGS, is_control=False, is_sv=False)

    preprocess.detect_sequence_error_reads(ARGS, is_control=False)
    preprocess.split_fastq_by_sequence_error(ARGS, is_control=False)
    preprocess.replace_midsv_without_sequence_errors(ARGS)

    preprocess.extract_knockin_loci(ARGS.tempdir, ARGS.sample_name)
    preprocess.cache_mutation_loci(ARGS, is_control=False)

    ###########################################################
    # Detect unexpected SV alleles on the control coordinate.
    ###########################################################
    preprocess.detect_sv_alleles(
        ARGS.tempdir, ARGS.sample_name, ARGS.control_name, ARGS.fasta_alleles, sv_type="insertion"
    )
    preprocess.detect_sv_alleles(
        ARGS.tempdir, ARGS.sample_name, ARGS.control_name, ARGS.fasta_alleles, sv_type="deletion"
    )
    preprocess.detect_sv_alleles(
        ARGS.tempdir, ARGS.sample_name, ARGS.control_name, ARGS.fasta_alleles, sv_type="inversion"
    )

    _add_predicted_sv_alleles(ARGS, original_alleles)

    ###########################################################
    # Score every read against every allele, then remap only the assigned reads.
    ###########################################################
    logger.info(f"Classify {arguments['sample']} with control-coordinate scoring...")
    allele_signatures = control_coordinate.build_allele_signatures(ARGS)
    control_midsv_records = control_coordinate.load_control_midsv_records(ARGS.tempdir, ARGS.sample_name)
    _classified_control, assignments = control_coordinate.classify_records_by_control_coordinate(
        control_midsv_records, allele_signatures, ARGS.no_filter
    )
    if not assignments:
        raise ValueError("No reads could be assigned by control-coordinate scoring.")
    fileio.save_pickle(
        assignments,
        Path(ARGS.tempdir, ARGS.sample_name, "classification", "control_coordinate_assignments.pickle"),
    )

    path_sample_fastq = Path(ARGS.tempdir, ARGS.sample_name, "fastq", f"{ARGS.sample_name}.fastq.gz")
    path_assigned_fastq = Path(ARGS.tempdir, ARGS.sample_name, "fastq", "control_coordinate_assignments")
    paths_fastq_by_allele = control_coordinate.split_fastq_by_assignment(
        path_sample_fastq, assignments, path_assigned_fastq
    )
    present_alleles = set(paths_fastq_by_allele)

    _copy_user_fastas_to_sample(ARGS, present_alleles)
    paths_sample_fasta = _paths_fasta_by_allele(ARGS.tempdir, ARGS.sample_name, present_alleles)
    control_coordinate.write_assigned_sam_files(
        ARGS, paths_fastq_by_allele, paths_sample_fasta, ARGS.sample_name, is_control_sv=False
    )
    preprocess.generate_midsv(ARGS, is_control=False, is_sv=False, target_alleles=present_alleles)

    path_control_fastq = Path(ARGS.tempdir, ARGS.control_name, "fastq", f"{ARGS.control_name}.fastq.gz")
    user_present_alleles = present_alleles & original_alleles
    sv_present_alleles = present_alleles - original_alleles

    if user_present_alleles:
        paths_control_fastq = dict.fromkeys(user_present_alleles, path_control_fastq)
        paths_control_fasta = _paths_fasta_by_allele(ARGS.tempdir, ARGS.control_name, user_present_alleles)
        control_coordinate.write_assigned_sam_files(
            ARGS, paths_control_fastq, paths_control_fasta, ARGS.control_name, is_control_sv=False
        )
        preprocess.generate_midsv(ARGS, is_control=True, is_sv=False, target_alleles=user_present_alleles)

    if sv_present_alleles:
        paths_control_fastq = dict.fromkeys(sv_present_alleles, path_control_fastq)
        paths_sv_fasta = _paths_fasta_by_allele(ARGS.tempdir, ARGS.sample_name, sv_present_alleles)
        control_coordinate.write_assigned_sam_files(
            ARGS, paths_control_fastq, paths_sv_fasta, ARGS.control_name, is_control_sv=True
        )
        preprocess.generate_midsv(ARGS, is_control=True, is_sv=True, target_alleles=sv_present_alleles)

    preprocess.cache_mutation_loci(ARGS, is_control=True)
    preprocess.extract_knockin_loci(ARGS.tempdir, ARGS.sample_name)
    preprocess.cache_mutation_loci(ARGS, is_control=False)

    allele_handler.save_allele_name_map(ARGS.tempdir, ARGS.fasta_alleles.keys())
    fileio.save_pickle(ARGS.fasta_alleles, Path(ARGS.tempdir, ARGS.sample_name, "fasta", "fasta_alleles.pickle"))

    classif_sample = control_coordinate.load_classified_midsv_for_assignments(ARGS, assignments)
    if not classif_sample:
        raise ValueError("No assigned reads produced allele-specific MIDSV records.")
    _complete_sample_analysis(arguments, ARGS, classif_sample)


def _add_predicted_sv_alleles(ARGS: FormattedInputs, original_alleles: set[str]) -> None:
    predefined_fasta = {
        str(Path(ARGS.tempdir, ARGS.sample_name, "fasta", f"{allele_handler.to_allele_key(allele)}.fasta"))
        for allele in original_alleles
    }
    paths_all_fasta = {str(path) for path in Path(ARGS.tempdir, ARGS.sample_name, "fasta").glob("*.fasta")}
    for path_fasta in paths_all_fasta - predefined_fasta:
        record = next(fileio.read_fasta(path_fasta))
        ARGS.fasta_alleles[record["identifier"]] = record["sequence"]


def _copy_user_fastas_to_sample(ARGS: FormattedInputs, alleles: set[str]) -> None:
    for allele in alleles:
        allele_key = allele_handler.to_allele_key(allele)
        path_source = Path(ARGS.tempdir, ARGS.control_name, "fasta", f"{allele_key}.fasta")
        path_destination = Path(ARGS.tempdir, ARGS.sample_name, "fasta", f"{allele_key}.fasta")
        if path_source.exists() and not path_destination.exists():
            shutil.copy(path_source, path_destination)


def _paths_fasta_by_allele(tempdir: Path, name: str, alleles: set[str]) -> dict[str, Path]:
    return {
        allele: Path(tempdir, name, "fasta", f"{allele_handler.to_allele_key(allele)}.fasta") for allele in alleles
    }


def _complete_sample_analysis(arguments: dict, ARGS: FormattedInputs, classif_sample: list[dict]) -> None:
    fileio.save_pickle(classif_sample, Path(ARGS.tempdir, ARGS.sample_name, "classification", "classif_sample.pickle"))

    ########################################################################
    # Clustering
    ########################################################################

    logger.info(f"Clustering {arguments['sample']}...")

    labels = clustering.extract_labels(classif_sample, ARGS.tempdir, ARGS.sample_name, ARGS.control_name)
    clust_sample = clustering.add_labels(classif_sample, labels)
    clust_sample = clustering.add_readnum(clust_sample)
    clust_sample = clustering.add_percent(clust_sample)
    clust_sample = clustering.update_labels(clust_sample)

    fileio.save_pickle(clust_sample, Path(ARGS.tempdir, ARGS.sample_name, "clustering", "clust_sample.pickle"))

    ########################################################################
    # Consensus call
    ########################################################################

    logger.info(f"Consensus calling of {arguments['sample']}...")

    # Remove minor alleles with fewer than 5 reads or less than 0.5%
    clust_sample_removed = consensus.remove_minor_alleles(clust_sample, ARGS.no_filter)

    # Adjust the percentage to 100%
    clust_sample_removed = consensus.scale_percentage(clust_sample_removed)

    # Downsampling to 1000 reads in each LABEL
    clust_downsampled = consensus.downsample_by_label(clust_sample_removed, 1000)

    consensus.cache_mutation_loci(ARGS, clust_downsampled, ARGS.no_filter)

    cons_percentages, cons_sequences, cons_midsv_tags, label_before_to_after = consensus.call_consensus(
        ARGS.tempdir, ARGS.sample_name, ARGS.fasta_alleles, clust_downsampled
    )

    map_label_name, map_name_allele = consensus.call_allele_name(cons_sequences, ARGS.fasta_alleles)
    cons_percentages = consensus.update_key_by_allele_name(cons_percentages, map_label_name)
    cons_sequences = consensus.update_key_by_allele_name(cons_sequences, map_label_name)
    cons_midsv_tags = consensus.update_key_by_allele_name(cons_midsv_tags, map_label_name)

    fileio.save_pickle(cons_percentages, Path(ARGS.tempdir, ARGS.sample_name, "consensus", "cons_percentages.pickle"))
    fileio.save_pickle(cons_sequences, Path(ARGS.tempdir, ARGS.sample_name, "consensus", "cons_sequences.pickle"))
    fileio.save_pickle(cons_midsv_tags, Path(ARGS.tempdir, ARGS.sample_name, "consensus", "cons_midsv_tags.pickle"))
    fileio.save_pickle(map_label_name, Path(ARGS.tempdir, ARGS.sample_name, "consensus", "map_label_name.pickle"))
    fileio.save_pickle(map_name_allele, Path(ARGS.tempdir, ARGS.sample_name, "consensus", "map_name_allele.pickle"))

    ########################################################################
    # Output Report：RESULT/FASTA/HTML/BAM
    ########################################################################

    logger.info(f"Output reports of {arguments['sample']}...")

    # RESULT
    RESULT_SAMPLE = consensus.update_label_percent_readnum_name(
        clust_sample_removed, map_label_name, label_before_to_after
    )
    RESULT_SAMPLE.sort(key=lambda x: x["LABEL"])

    fileio.write_jsonl(RESULT_SAMPLE, Path(ARGS.tempdir, "result", f"{ARGS.sample_name}.jsonl"))

    # FASTA
    report.sequence_exporter.export_to_fasta(ARGS.tempdir, ARGS.sample_name, cons_sequences)
    report.sequence_exporter.export_reference_to_fasta(ARGS.tempdir, ARGS.sample_name)

    # HTML
    report.sequence_exporter.export_to_html(
        ARGS.tempdir, ARGS.sample_name, ARGS.fasta_alleles, cons_midsv_tags, map_name_allele
    )

    # VCF
    report.vcf_exporter.export_to_vcf(ARGS.tempdir, ARGS.sample_name, ARGS.genome_coordinates, cons_midsv_tags)

    # BAM
    report.bam_exporter.export_to_bam(
        ARGS.tempdir, ARGS.sample_name, ARGS.genome_coordinates, ARGS.threads, RESULT_SAMPLE
    )
    for path_bam_igvjs in Path(ARGS.tempdir, "cache", ".igvjs").glob(f"{ARGS.control_name}_control.bam*"):
        shutil.copy(path_bam_igvjs, Path(ARGS.tempdir, "report", ".igvjs", ARGS.sample_name))

    ###########################################################
    # Finish call
    ###########################################################

    logger.info(f"\N{TEACUP WITHOUT HANDLE} {arguments['sample']} is finished!")
