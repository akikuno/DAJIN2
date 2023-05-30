from __future__ import annotations

import sys

from DAJIN2.core import core_execute
from DAJIN2.postprocess import report
from DAJIN2.preprocess.validate_inputs import (validate_files,
                                               validate_genome_and_fetch_urls)


def single_execute(arguments: dict[str]):
    ################################################################################
    # Validate contents in the batch file
    ################################################################################
    validate_files(args["sample"], args["control"], args["allele"])
    if "genome" in args:
        UCSC_URL, GOLDENPATH_URL = validate_genome_and_fetch_urls(args["genome"])
    else:
        UCSC_URL, GOLDENPATH_URL = None, None

    print(f"{arguments['control']} is now processing...", file=sys.stderr)
    core_execute.execute_control(arguments)
    print(f"{arguments['sample']} is now processing...", file=sys.stderr)
    core_execute.execute_sample(arguments)
    name = arguments["name"]
    report.report(name)
    print(f"Finished! Open DAJINResults/{name} to see the report.", file=sys.stderr)
