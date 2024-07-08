from __future__ import annotations

import datetime
import logging
import sys
import warnings
from pathlib import Path

from sklearn.exceptions import ConvergenceWarning

DAJIN_VERSION = "0.5.2"
DAJIN_RESULTS_DIR = Path("DAJIN_Results")
TEMP_ROOT_DIR = Path(DAJIN_RESULTS_DIR, ".tempdir")


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
    return Path(f"DAJIN2_log_{current_time}.txt")


def set_logging(path_logfile: Path) -> logging.Logger:
    log_format = "%(asctime)s, %(levelname)s, %(message)s"
    date_format = "%Y-%m-%d %H:%M:%S"

    stderr_handler = logging.StreamHandler()
    stderr_handler.setFormatter(logging.Formatter(log_format, datefmt=date_format))

    file_handler = DeferredFileHandler(path_logfile)
    file_handler.setFormatter(logging.Formatter(log_format, datefmt=date_format))

    logging.basicConfig(
        level=logging.INFO,
        handlers=[stderr_handler, file_handler],
    )

    # log uncaught exceptions
    def handle_uncaught_exception(exc_type, exc_value, exc_traceback):
        logging.error("Catch an Exception. Traceback:", exc_info=(exc_type, exc_value, exc_traceback))

    sys.excepthook = handle_uncaught_exception

    return logging.getLogger(__name__)


def reset_logging() -> None:
    """
    Reset the logging system to its default state.
    """
    # Remove all existing handlers
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)


def set_warnings_ignore() -> None:
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", category=UserWarning)
    warnings.filterwarnings("ignore", category=ConvergenceWarning)
