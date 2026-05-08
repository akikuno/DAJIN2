from __future__ import annotations

from typing import Any

from .._lazy_import import LazyExportMap, list_lazy_exports, load_lazy_export

_LAZY_EXPORTS: LazyExportMap = {
    "add_labels": (".appender", "add_labels"),
    "add_percent": (".appender", "add_percent"),
    "add_readnum": (".appender", "add_readnum"),
    "extract_labels": (".label_extractor", "extract_labels"),
    "update_labels": (".label_updator", "update_labels"),
}

__all__ = sorted(_LAZY_EXPORTS)


def __getattr__(name: str) -> Any:
    return load_lazy_export(globals(), _LAZY_EXPORTS, name)


def __dir__() -> list[str]:
    return list_lazy_exports(globals(), _LAZY_EXPORTS)
