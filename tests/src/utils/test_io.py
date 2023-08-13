from src.DAJIN2.utils import io
import os
import json

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
    return test == answer


########################################################################
# Input/Output Json
########################################################################


def test_write_jsonl():
    data_list = [{"name": "Alice", "age": 30}, {"name": "Bob", "age": 25}]
    test_filename = "test_output.json"
    io.write_jsonl(test_filename, data_list)
    # Verify if the file has been written correctly
    with open(test_filename, "r") as f:
        lines = f.readlines()
        assert len(lines) == 2
        loaded_data1 = json.loads(lines[0].strip())
        loaded_data2 = json.loads(lines[1].strip())
        assert loaded_data1 == {"name": "Alice", "age": 30}
        assert loaded_data2 == {"name": "Bob", "age": 25}
    # remove test file
    os.remove(test_filename)
