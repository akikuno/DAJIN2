from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from src.DAJIN2.utils import io


def test_sanitize_name_with_valid_path():
    assert io.sanitize_name("valid_name") == "valid_name"
    assert io.sanitize_name(Path("valid/name")) == "valid-name"


def test_sanitize_name_with_invalid_characters():
    assert io.sanitize_name("inva/lid:name?") == "inva-lid-name-"
    assert io.sanitize_name(Path("/inva>lid|name.")) == "-inva-lid-name-"


def test_sanitize_name_with_whitespace():
    assert io.sanitize_name("  leading space") == "leading-space"
    assert io.sanitize_name("trailing space ") == "trailing-space"


def test_sanitize_name_with_empty_string():
    with pytest.raises(ValueError) as e:
        io.sanitize_name(" ")
    assert str(e.value) == "Provided name is empty or consists only of whitespace"


def test_sanitize_name_with_empty_path():
    with pytest.raises(ValueError) as e:
        io.sanitize_name("")
    assert str(e.value) == "Provided name is empty or consists only of whitespace"


########################################################################
# Convert Path
########################################################################


def test_convert_to_posix_winpath():
    path = r"C:\Windows\System32"
    test = io.convert_to_posix(path)
    answer = "/mnt/c/Windows/System32"
    assert test == answer


def test_convert_to_posix_posixpath():
    path = r"/mnt/c/Windows/System32"
    test = io.convert_to_posix(path)
    answer = "/mnt/c/Windows/System32"
    assert test == answer


def test_convert_to_posix():
    path = "example\\Desktop\\test.txt"
    test = io.convert_to_posix(path)
    answer = "example/Desktop/test.txt"
    assert test == answer


########################################################################
# Input/Output Json
########################################################################


def test_write_jsonl():
    data_list = [{"name": "Alice", "age": 30}, {"name": "Bob", "age": 25}]
    test_filename = "test_output.json"
    io.write_jsonl(file_path=test_filename, data=data_list)
    # Verify if the file has been written correctly
    with open(test_filename) as f:
        lines = f.readlines()
        assert len(lines) == 2
        loaded_data1 = json.loads(lines[0].strip())
        loaded_data2 = json.loads(lines[1].strip())
        assert loaded_data1 == {"name": "Alice", "age": 30}
        assert loaded_data2 == {"name": "Bob", "age": 25}
    # remove test file
    os.remove(test_filename)


########################################################################
# is_gip_file
########################################################################


@pytest.mark.parametrize(
    "path_file, expected",
    [
        ("tests/data/utils/io/test.fasta.gz", True),
        ("tests/data/utils/io/test.fasta", False),
    ],
)
def test_is_gzip_file(path_file, expected):
    assert io.is_gzip_file(path_file) == expected


# Utility function to create a temporary file and write some data into it
def create_temp_file(tmp_path, filename, content):
    file_path = tmp_path / filename
    with open(file_path, "wb") as f:
        f.write(content)
    return file_path


# Test cases
def test_count_newlines_empty_file(tmp_path):
    file_path = create_temp_file(tmp_path, "empty.txt", b"")
    assert io.count_newlines(file_path) == 0


def test_count_newlines_single_line_no_newline(tmp_path):
    file_path = create_temp_file(tmp_path, "single_line.txt", b"Hello, world!\n")
    assert io.count_newlines(file_path) == 1


def test_count_newlines_multiple_lines(tmp_path):
    file_path = create_temp_file(tmp_path, "multiple_lines.txt", b"Hello, world!\nHow are you?\nI am good.\n")
    assert io.count_newlines(file_path) == 3
