from __future__ import annotations

from .core import core
from .postprocess import report


def execute(arguments: dict[str]):
    core.execute_preprocess(arguments)
    core.execute_control(arguments)
    core.execute_sample(arguments)
    name = arguments["name"]
    report.report(name)
