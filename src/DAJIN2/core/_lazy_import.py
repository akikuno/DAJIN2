from __future__ import annotations

from importlib import import_module
from typing import Any, TypeAlias

LazyExportMap: TypeAlias = dict[str, tuple[str, str | None]]


def load_lazy_export(module_globals: dict[str, Any], exports: LazyExportMap, name: str) -> Any:
    try:
        module_path, attribute_name = exports[name]
    except KeyError:
        module_name = module_globals.get("__name__", "module")
        raise AttributeError(f"module {module_name!r} has no attribute {name!r}") from None

    module = import_module(module_path, package=module_globals["__name__"])
    value = module if attribute_name is None else getattr(module, attribute_name)
    module_globals[name] = value
    return value


def list_lazy_exports(module_globals: dict[str, Any], exports: LazyExportMap) -> list[str]:
    return sorted(set(module_globals) | set(exports))
