from __future__ import annotations

import sys
import traceback
import logging

from typing import Generator
from itertools import islice
from multiprocessing import Process, Queue


def _batched(iterable, chunk_size: int) -> Generator(tuple):
    iterator = iter(iterable)
    while True:
        chunk = tuple(islice(iterator, chunk_size))
        if not chunk:
            return
        yield chunk


def run(function, arguments: list[dict], num_workers: int = 1) -> None:
    logger = logging.getLogger(__name__)

    q = Queue()

    def target(arg):
        try:
            function(arg)
        except Exception:
            q.put(traceback.format_exc())
            sys.exit(1)

    arguments_batched = _batched(arguments, num_workers)
    for args in arguments_batched:
        processes = [Process(target=target, args=(arg,)) for arg in args]

        for p in processes:
            p.start()

        for p in processes:
            p.join()
            if p.exitcode == 1:
                logger.error(f"An unexpected error occurred: {q.get()}")
                sys.exit(1)
