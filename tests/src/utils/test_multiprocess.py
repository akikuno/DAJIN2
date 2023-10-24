from __future__ import annotations
import pytest

import os
import re
import logging
import tempfile
from src.DAJIN2.utils import multiprocess


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
# run_multiprocess
###########################################################


def _setup_logging(log_file_path):
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    logging.basicConfig(filename=log_file_path, level=logging.INFO)


def _dummy_function(args: dict[str]):
    n = args["value"]
    temp_file_path = args["path"]
    _setup_logging("/tmp/multiprocess.log")  # 各プロセスでのログ設定
    with open(temp_file_path, "a") as f:
        f.write(str(n) + "\n")


@pytest.mark.slow
def test_run_multiprocess():
    # Use tempfile to create a temporary file
    with tempfile.NamedTemporaryFile(delete=True) as temp_file:
        temp_file_path = temp_file.name

        # Arguments to run the _dummy_function
        arguments = [
            {"value": 1, "path": temp_file_path},
            {"value": 2, "path": temp_file_path},
            {"value": 3, "path": temp_file_path},
            {"value": 4, "path": temp_file_path},
            {"value": 5, "path": temp_file_path},
            {"value": 6, "path": temp_file_path},
            {"value": 7, "path": temp_file_path},
            {"value": 8, "path": temp_file_path},
            {"value": 9, "path": temp_file_path},
            {"value": 10, "path": temp_file_path},
        ]

        # Use the run_multiprocess function
        multiprocess.run(_dummy_function, arguments, num_workers=3)

        # Check if the _dummy_function wrote to the file correctly
        with open(temp_file_path, "r") as f:
            lines = f.readlines()

        assert set(map(int, map(str.strip, lines))) == {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}


###########################################################
# test logging of run_multiprocess
###########################################################


def _dummy_function_logging(args: dict[str]):
    n = args["value"]
    temp_file_path = args["path"]
    _setup_logging("/tmp/multiprocess.log")  # 各プロセスでのログ設定
    with open(temp_file_path, "a") as f:
        f.write(str(n) + "\n")


def _dummy_function_that_fails(arg):
    raise ValueError("This is a simulated error!")


@pytest.mark.slow
def test_logging_of_run_multiprocess():
    # ログの設定を一時的に変更
    with tempfile.NamedTemporaryFile(delete=True, mode="w+") as temp_log_file:
        log_file_path = temp_log_file.name
        _setup_logging(log_file_path)  # ログ設定のリセットと再設定

        # run_multiprocess関数を実行
        arguments = [{"value": 1, "path": "/tmp/dummy_path"}, {"value": 2, "path": "/tmp/dummy_path"}]
        multiprocess.run(_dummy_function_logging, arguments, num_workers=2)

        # ログファイルから内容を読み込む
        temp_log_file.seek(0)
        logs = temp_log_file.read()

        # 期待するログメッセージが存在するかどうかを確認
        assert "Starting process Process-" in logs


@pytest.mark.slow
def test_run_multiprocess_failure():
    arguments = [{"dummy_arg": 1}, {"dummy_arg": 2}]

    # run_multiprocess関数がSystemExitを引き起こすことを確認
    with pytest.raises(SystemExit) as exc_info:
        multiprocess.run(_dummy_function_that_fails, arguments, num_workers=2)

    # exitcode 1が期待される
    assert exc_info.value.code == 1


###########################################################
# remove log files
###########################################################
def delete_log_files(directory="."):
    # パターンに一致するファイル名を探す正規表現
    pattern = r"^\d{14}_DAJIN2\.log$"
    # 指定されたディレクトリ内のファイルをリストアップ
    for filename in os.listdir(directory):
        if re.match(pattern, filename):
            filepath = os.path.join(directory, filename)
            try:
                os.remove(filepath)
                print(f"Deleted: {filepath}")
            except Exception as e:
                print(f"Error deleting {filepath}: {e}")


delete_log_files()
