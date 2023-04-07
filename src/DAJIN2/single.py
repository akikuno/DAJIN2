from __future__ import annotations

from DAJIN2.core import core_execute
from DAJIN2.postprocess import report


def single_execute(arguments: dict[str]):
    core_execute.execute_control(arguments)
    core_execute.execute_sample(arguments)
    name = arguments["name"]
    report.report(name)
