from __future__ import annotations

from typing import Any

from .._lazy_import import LazyExportMap, list_lazy_exports, load_lazy_export

_LAZY_EXPORTS: LazyExportMap = {
    "cache_mutation_loci": (".consensus_mutation_analyzer", "cache_mutation_loci"),
    "call_allele_name": (".consensus_formatter", "call_allele_name"),
    "call_consensus": (".consensus", "call_consensus"),
    "downsample_by_label": (".clust_formatter", "downsample_by_label"),
    "remove_minor_alleles": (".clust_formatter", "remove_minor_alleles"),
    "scale_percentage": (".consensus_formatter", "scale_percentage"),
    "update_key_by_allele_name": (".consensus_formatter", "update_key_by_allele_name"),
    "update_label_percent_readnum_name": (".consensus_formatter", "update_label_percent_readnum_name"),
}

__all__ = sorted(_LAZY_EXPORTS)


def __getattr__(name: str) -> Any:
    return load_lazy_export(globals(), _LAZY_EXPORTS, name)


def __dir__() -> list[str]:
    return list_lazy_exports(globals(), _LAZY_EXPORTS)
