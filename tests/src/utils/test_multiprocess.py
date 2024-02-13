from __future__ import annotations

import pytest

import tempfile
from src.DAJIN2.utils import multiprocess

from multiprocessing import Queue

###########################################################
# generate_chunks
###########################################################


def test_generate_chunks():
    # Basic test to check if it chunks correctly
    iterable = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    chunk_size = 3
    result = list(multiprocess.generate_chunks(iterable, chunk_size))
    assert result == [(1, 2, 3), (4, 5, 6), (7, 8, 9), (10,)]


def test_generate_chunks_larger():
    # Chunk size larger than iterable
    iterable = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    chunk_size = 20
    result = list(multiprocess.generate_chunks(iterable, chunk_size))
    assert result == [(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)]


def test_generate_chunks_empty():
    # Empty iterable
    iterable = []
    chunk_size = 3
    result = list(multiprocess.generate_chunks(iterable, chunk_size))
    assert result == []


def test_generate_chunks_one():
    # Chunk size of 1
    iterable = [1, 2, 3, 4, 5]
    chunk_size = 1
    result = list(multiprocess.generate_chunks(iterable, chunk_size))
    assert result == [(1,), (2,), (3,), (4,), (5,)]


###########################################################
# target
###########################################################


def sample_function(arg: dict) -> None:
    """Sample function to demonstrate how the `run` function works."""
    if "error" in arg["sample"]:
        raise ValueError("Sample error!")


def test_target_without_exception():
    queue = Queue()
    multiprocess.target(sample_function, {"sample": "normal"}, queue)
    assert queue.empty()  # Queue should be empty if no error


def test_target_with_exception():
    queue = Queue()
    with pytest.raises(ValueError):
        multiprocess.target(sample_function, {"sample": "error"}, queue)

    error_message = queue.get()
    assert "Sample error!" in error_message
    assert "An unexpected error occurred at error" in error_message


###########################################################
# run_multiprocess
###########################################################


def test_run_without_exception(caplog):
    arguments = [{"sample": "normal1"}, {"sample": "normal2"}]
    multiprocess.run(sample_function, arguments, num_workers=2)
    assert "An unexpected error occurred" not in caplog.text


def test_run_with_exception(caplog):
    arguments = [{"sample": "normal1"}, {"sample": "error"}]
    multiprocess.run(sample_function, arguments, num_workers=2)
    assert "Sample error!" in caplog.text


def test_run_with_multiple_exceptions(caplog):
    arguments = [{"sample": "error1"}, {"sample": "error2"}]
    multiprocess.run(sample_function, arguments, num_workers=2)
    assert "Sample error!" in caplog.text
    assert (
        "An unexpected error occurred at error1" in caplog.text
        or "An unexpected error occurred at error2" in caplog.text
    )


def write_value_to_file(args: dict) -> None:
    """
    Writes the given value to the specified file and sets up logging.
    """
    value = args["value"]
    file_path = args["path"]

    with open(file_path, "a") as f:
        f.write(str(value) + "\n")


@pytest.mark.slow
def test_multiprocessing_execution():
    """
    Test if the function `write_value_to_file` can be executed in parallel using multiple processes.
    """
    # Create a temporary file to store the values
    with tempfile.NamedTemporaryFile(delete=True) as temp_file:
        file_path = temp_file.name

        # Prepare arguments for parallel execution
        values_to_write = [{"value": i, "path": file_path, "sample": "test"} for i in range(1, 11)]

        # Execute the function in parallel using 3 worker processes
        multiprocess.run(write_value_to_file, values_to_write, num_workers=3)

        # Verify if all values are written correctly
        with open(file_path, "r") as f:
            written_values = set(map(int, map(str.strip, f.readlines())))

        assert written_values == set(range(1, 11))
