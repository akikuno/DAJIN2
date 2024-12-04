from __future__ import annotations

import io
import logging
import sys
import traceback
from collections.abc import Callable, Iterator
from itertools import islice
from multiprocessing import Process, Queue


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


def handle_exception_and_enqueue(queue: Queue[str], arg: dict) -> None:
    """
    Handle the exception by placing an error message in the queue.
    """
    error_message = f"{get_error_message_prefix(arg)}\n{traceback.format_exc()}"
    queue.put(error_message)


def target(function: Callable, arg: dict, queue: Queue[str]) -> None:
    """Run a function with a single argument and handle any exceptions."""
    try:
        suppress_stderr(function, arg)
    except Exception:
        handle_exception_and_enqueue(queue, arg)
        raise


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
                name=get_error_message_prefix(arg),
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
