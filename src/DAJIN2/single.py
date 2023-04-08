from __future__ import annotations

import sys
from DAJIN2.core import core_execute
from DAJIN2.postprocess import report


def single_execute(arguments: dict[str]):
    print(f"{arguments['control']} is now processing...", file=sys.stderr)
    core_execute.execute_control(arguments)
    print(f"{arguments['sample']} is now processing...", file=sys.stderr)
    core_execute.execute_sample(arguments)
    name = arguments["name"]
    report.report(name)
    print(f"Finished! Open DAJINResults/{name} to see the report.", file=sys.stderr)
