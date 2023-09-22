from __future__ import annotations

import io
import sys
import logging
import traceback

from itertools import islice
from multiprocessing import Process, Queue

from typing import Generator
from collections.abc import Callable


def generate_chunks(iterable, chunk_size: int = 1) -> Generator[tuple]:
    """Yield successive n-sized chunks from iterable."""
    iterator = iter(iterable)
    while True:
        chunk = tuple(islice(iterator, chunk_size))
        if not chunk:
            return
        yield chunk


def target(function: Callable, arg: dict, queue: Queue[str]) -> None:
    """Run a function with a single argument, and handle any exceptions."""
    try:
        # Redirect stderr to suppress traceback printing
        original_stderr = sys.stderr
        sys.stderr = io.StringIO()
        function(arg)
    except Exception:
        # Reset stderr back to original
        sys.stderr = original_stderr
        queue.put(f"An unexpected error occurred at {arg['sample']}\n{traceback.format_exc()}")
        raise  # Re-raise the exception to be caught in the parent process
    finally:
        # Ensure stderr is reset back to original even if no error occurs
        sys.stderr = original_stderr


def run(function: Callable, arguments: list[dict], num_workers: int = 1) -> None:
    """Run a function in parallel over a list of arguments."""
    logger = logging.getLogger(__name__)

    queue = Queue()

    arguments_chunked = generate_chunks(arguments, num_workers)
    for args in arguments_chunked:
        processes = [
            Process(
                target=target,
                args=(function, arg, queue),
                name=f"An unexpected error occurred at {arg['sample']}",
            )
            for arg in args
        ]

        for process in processes:
            process.start()

        for process in processes:
            process.join()
            if process.exitcode == 1:
                logger.error(queue.get())
                break
