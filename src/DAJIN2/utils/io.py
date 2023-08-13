from __future__ import annotations
import os
import wslPath
import json
import pandas as pd
import csv
from typing import Generator


###########################################################
# convert path to POSIX format
###########################################################


def convert_to_posix(path: str) -> str:
    if wslPath.is_windows_path(path):
        path = wslPath.to_posix(path)
    return path


###########################################################
# Load batch file
###########################################################


def check_excel_or_csv(filepath: str) -> str:
    """Check if the file is an Excel or CSV file. Raise error for other types."""
    _, ext = os.path.splitext(filepath)
    if ext in [".xlsx", ".xls"]:
        return "excel"
    elif ext == ".csv":
        return "csv"
    else:
        raise ValueError("The provided file must be either an Excel or CSV file.")


def load_from_excel(path: str) -> list:
    """Load data from an Excel file and return as a list."""
    df = pd.read_excel(path)
    return [df.columns.to_list()] + df.values.tolist()


def load_from_csv(path: str) -> list:
    """Load data from a CSV file and return as a list."""
    with open(path, "r") as f:
        return [row for row in csv.reader(f, skipinitialspace=True, delimiter=",")]


def load_file(path_batchfile: str) -> list:
    """Load data from either an Excel or CSV file."""
    file_type = check_excel_or_csv(path_batchfile)
    if file_type == "excel":
        return load_from_excel(path_batchfile)
    elif file_type == "csv":
        return load_from_csv(path_batchfile)


###########################################################
# Input/Output Json
###########################################################


def read_jsonl(path: str) -> Generator:
    with open(path, "r") as f:
        for line in f:
            yield json.loads(line)


def write_jsonl(path: str, data: Generator) -> None:
    with open(path, "w") as f:
        for d in data:
            json_str = json.dumps(d)
            f.write(json_str + "\n")
