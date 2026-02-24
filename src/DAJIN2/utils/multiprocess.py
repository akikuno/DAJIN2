from __future__ import annotations

import io
import logging
import os
import queue
import sys
import time
import traceback
from collections.abc import Callable, Iterator
from itertools import islice
from multiprocessing import Queue, get_context
from pathlib import Path

from DAJIN2.utils import config


def get_error_message_prefix(arg: dict) -> str:
    """Generate the error message prefix based on the argument's sample value."""
    return f"An unexpected error occurred at {arg['sample']}"


def generate_chunks(iterable, chunk_size: int = 1) -> Iterator[tuple]:
    """Yield successive n-sized chunks from iterable."""
    iterator = iter(iterable)
    while True:
        chunk = tuple(islice(iterator, chunk_size))
        if not chunk:
            return
        yield chunk


def suppress_stderr(function: Callable, *args, **kwargs) -> None:
    """
    Suppress stderr and execute the function.
    If an error occurs, revert the stderr and raise the exception.
    """
    original_stderr = sys.stderr
    sys.stderr = io.StringIO()

    try:
        function(*args, **kwargs)
    finally:
        sys.stderr = original_stderr


class ProgressQueueLogHandler(logging.Handler):
    """Forward child process logs to the GUI progress queue."""

    def __init__(self, progress_log_queue: Queue[str]):
        super().__init__()
        self.progress_log_queue = progress_log_queue

    def emit(self, record: logging.LogRecord) -> None:
        try:
            self.progress_log_queue.put_nowait(self.format(record))
        except Exception:
            # Never stop analysis because of logging relay failures.
            pass


def setup_child_logging(progress_log_queue: Queue[str] | None = None) -> None:
    root_logger = logging.getLogger()
    level_name = os.environ.get("DAJIN2_LOGLEVEL", "INFO").upper()
    level = getattr(logging, level_name, logging.INFO)

    if not root_logger.handlers:
        path_logfile = os.environ.get("DAJIN2_LOGFILE")
        if path_logfile:
            config.set_logging(Path(path_logfile), level=level)
        else:
            logging.basicConfig(level=level)

    if progress_log_queue is None:
        return

    if any(isinstance(handler, ProgressQueueLogHandler) for handler in root_logger.handlers):
        return

    progress_handler = ProgressQueueLogHandler(progress_log_queue)
    progress_handler.setLevel(level)
    formatter = logging.Formatter("%(asctime)s, %(levelname)s, %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    progress_handler.setFormatter(formatter)
    root_logger.addHandler(progress_handler)


def handle_exception_and_enqueue(queue: Queue[str], arg: dict) -> None:
    """
    Handle the exception by placing an error message in the queue.
    """
    error_message = f"{get_error_message_prefix(arg)}\n{traceback.format_exc()}"
    queue.put(error_message)


def target(function: Callable, arg: dict, queue: Queue[str], progress_log_queue: Queue[str] | None = None) -> None:
    """Run a function with a single argument and handle any exceptions."""
    setup_child_logging(progress_log_queue)
    try:
        suppress_stderr(function, arg)
    except Exception:
        handle_exception_and_enqueue(queue, arg)
        raise


def enqueue_progress_logs(log_queue: Queue[str] | None, progress_queue) -> None:
    if log_queue is None or progress_queue is None:
        return

    while True:
        try:
            message = log_queue.get_nowait()
            progress_queue.put({"status": "log", "message": message, "timestamp": time.time()})
        except queue.Empty:
            break


def run(function: Callable, arguments: list[dict], num_workers: int = 1, progress_queue=None) -> None:
    """Run a function in parallel over a list of arguments."""
    logger = logging.getLogger(__name__)

    ctx = get_context("spawn")
    error_queue = ctx.Queue()
    log_queue = ctx.Queue() if progress_queue is not None else None

    arguments_chunked = generate_chunks(arguments, num_workers)
    for args in arguments_chunked:
        processes = [
            ctx.Process(
                target=target,
                args=(function, arg, error_queue, log_queue),
                name=get_error_message_prefix(arg),
            )
            for arg in args
        ]

        for process in processes:
            process.start()

        while any(process.is_alive() for process in processes):
            enqueue_progress_logs(log_queue, progress_queue)
            for process in processes:
                process.join(timeout=0.1)

        enqueue_progress_logs(log_queue, progress_queue)

        for process in processes:
            process.join()
            if process.exitcode == 1:
                logger.error(error_queue.get())
                break
