from __future__ import annotations

import os
import logging
import warnings
import datetime

from pathlib import Path

from sklearn.exceptions import ConvergenceWarning


DAJIN_RESULTS_DIR = Path("DAJIN_Results")
TEMP_ROOT_DIR = Path(DAJIN_RESULTS_DIR, ".tempdir")


def set_single_threaded_blas():
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"


def set_logging():
    current_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

    stderr_handler = logging.StreamHandler()
    stderr_handler.setFormatter(logging.Formatter("%(message)s"))

    file_handler = logging.FileHandler(f"{current_time}_DAJIN2.log")
    file_handler.setFormatter(logging.Formatter("%(message)s"))

    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] [%(levelname)s] [%(processName)s] %(message)s",
        handlers=[stderr_handler, file_handler],
    )


def set_warnings():
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", category=ConvergenceWarning)
