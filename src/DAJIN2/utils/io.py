import gzip
import json


def fread(filepath: str) -> list:
    """
    read gzip or text file
    """
    with gzip.open(filepath, "rt") as f:
        try:
            out = f.read().splitlines()
        except gzip.BadGzipFile:
            with open(filepath) as f:
                out = f.read().splitlines()
    return out


def fwrite(infile: list, outfilepath: str) -> None:
    """
    save gzipped file
    """
    with gzip.open(outfilepath, "wt") as f:
        for i in infile:
            f.write("\n".join(i) + "\n")


def read_jsonl(filepath: str) -> list[dict]:
    """Read a JSONL file

    Args:
        filepath (str): Path to a file

    Returns:
        list[dict]: A list of dictionaries
    """
    dicts = []
    with open(filepath, "r") as f:
        for line in f:
            dicts.append(json.JSONDecoder().decode(line))
    return dicts


def write_jsonl(dicts: list[dict], filepath: str):
    """Write a dictionary to a JSONL

    Args:
        dicts (list[dict]): A dictionary to write to disk
        filepath (str): Path to a file to write to
    """
    with open(filepath, "w") as output:
        for d in dicts:
            json.dump(d, output)
            output.write("\n")
