from __future__ import annotations

from typing import Any

from .._lazy_import import LazyExportMap, list_lazy_exports, load_lazy_export

_LAZY_EXPORTS: LazyExportMap = {
    "bam_exporter": (".bam_exporter", None),
    "sequence_exporter": (".sequence_exporter", None),
    "vcf_exporter": (".vcf_exporter", None),
}

__all__ = sorted(_LAZY_EXPORTS)


def __getattr__(name: str) -> Any:
    return load_lazy_export(globals(), _LAZY_EXPORTS, name)


def __dir__() -> list[str]:
    return list_lazy_exports(globals(), _LAZY_EXPORTS)
