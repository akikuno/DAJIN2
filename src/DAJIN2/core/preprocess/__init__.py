from __future__ import annotations

from typing import Any

from .._lazy_import import LazyExportMap, list_lazy_exports, load_lazy_export

_LAZY_EXPORTS: LazyExportMap = {
    "cache_mutation_loci": (".mutation_processing.mutation_extractor", "cache_mutation_loci"),
    "create_report_directories": (".infrastructure.directory_manager", "create_report_directories"),
    "create_temporal_directories": (".infrastructure.directory_manager", "create_temporal_directories"),
    "detect_sequence_error_reads": (".error_correction.sequence_error_handler", "detect_sequence_error_reads"),
    "detect_sv_alleles": (".structural_variants.sv_detector", "detect_sv_alleles"),
    "extract_knockin_loci": (".mutation_processing.knockin_handler", "extract_knockin_loci"),
    "extract_sequence_errors_in_strand_biased_loci": (
        ".error_correction.strand_bias_handler",
        "extract_sequence_errors_in_strand_biased_loci",
    ),
    "fetch_chromosome_size": (".genome_coordinate.genome_fetcher", "fetch_chromosome_size"),
    "fetch_coordinates": (".genome_coordinate.genome_fetcher", "fetch_coordinates"),
    "format_inputs": (".infrastructure.input_formatter", "format_inputs"),
    "generate_midsv": (".alignment.midsv_caller", "generate_midsv"),
    "generate_sam": (".alignment.mapping", "generate_sam"),
    "replace_midsv_without_sequence_errors": (
        ".error_correction.sequence_error_handler",
        "replace_midsv_without_sequence_errors",
    ),
    "split_fastq_by_sequence_error": (".error_correction.sequence_error_handler", "split_fastq_by_sequence_error"),
}

__all__ = sorted(_LAZY_EXPORTS)


def __getattr__(name: str) -> Any:
    return load_lazy_export(globals(), _LAZY_EXPORTS, name)


def __dir__() -> list[str]:
    return list_lazy_exports(globals(), _LAZY_EXPORTS)
