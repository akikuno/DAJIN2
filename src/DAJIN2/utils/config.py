from __future__ import annotations

import os
import sys
import logging
import warnings
import datetime

from pathlib import Path

from sklearn.exceptions import ConvergenceWarning


DAJIN_RESULTS_DIR = Path("DAJIN_Results")
TEMP_ROOT_DIR = Path(DAJIN_RESULTS_DIR, ".tempdir")


def set_single_threaded_blas() -> None:
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"


class DeferredFileHandler(logging.FileHandler):
    def __init__(self, filename, mode="a", encoding=None, delay=True):
        # Setting delay to True to defer the file opening
        super().__init__(filename, mode, encoding, delay)

    def emit(self, record):
        # The file is actually opened when a log is emitted
        if self.stream is None:
            self.stream = self._open()
        super().emit(record)


def get_logfile() -> Path:
    current_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    return Path(f"{current_time}_DAJIN2.log")


def set_logging(path_logfile: Path) -> None:
    format = "%(asctime)s, %(levelname)s, %(message)s"
    datefmt = "%Y-%m-%d %H:%M:%S"

    stderr_handler = logging.StreamHandler()
    stderr_handler.setFormatter(logging.Formatter(format, datefmt=datefmt))

    file_handler = DeferredFileHandler(path_logfile)
    file_handler.setFormatter(logging.Formatter(format, datefmt=datefmt))

    logging.basicConfig(
        level=logging.INFO,
        handlers=[stderr_handler, file_handler],
    )

    # log uncaught exceptions
    def handle_uncaught_exception(exc_type, exc_value, exc_traceback):
        logging.error("Catch an Exception. Traceback:", exc_info=(exc_type, exc_value, exc_traceback))

    sys.excepthook = handle_uncaught_exception


def set_warnings() -> None:
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", category=ConvergenceWarning)
