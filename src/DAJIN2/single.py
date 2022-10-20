from __future__ import annotations
from .core import core
from .postprocess import report


def execute(arguments: dict[str]):
    core.execute(arguments)
    name = arguments["name"]
    report.report(name)
