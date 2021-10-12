import os
from pprint import pprint
from src.DAJIN2.preprocess import mids

test_dir = os.path.join(
    "tests", "data", "preprocess",
    "mids", "samfile_to_mids"
)

input_dir = os.path.join(test_dir, "input")
expectation_dir = os.path.join(test_dir, "expectation")

input_files = [os.path.join(input_dir, _) for _ in os.listdir(input_dir)]
expectation_files = [os.path.join(expectation_dir, _)
                     for _ in os.listdir(expectation_dir)]

pprint(input_files)
pprint(expectation_files)


def test_ins():
    infile = [s for s in input_files if "ins" in s][0]
    expfile = [s for s in expectation_files if "ins" in s][0]
    with open(expfile) as f:
        expected = f.readlines()
    assert mids.samfile_to_mids(infile) == expected


def test_del():
    infile = [s for s in input_files if "del.sam" in s][0]
    expfile = [s for s in expectation_files if "del.csv" in s][0]
    with open(expfile) as f:
        expected = f.readlines()
    assert mids.samfile_to_mids(infile) == expected


def test_sub():
    infile = [s for s in input_files if "sub.sam" in s][0]
    expfile = [s for s in expectation_files if "sub.csv" in s][0]
    with open(expfile) as f:
        expected = f.readlines()
    assert mids.samfile_to_mids(infile) == expected


def test_del_long():
    infile = [s for s in input_files if "del_long.sam" in s][0]
    expfile = [s for s in expectation_files if "del_long.csv" in s][0]
    with open(expfile) as f:
        expected = f.readlines()
    assert mids.samfile_to_mids(infile) == expected


def test_inv():
    infile = [s for s in input_files if "inv.sam" in s][0]
    expfile = [s for s in expectation_files if "inv.csv" in s][0]
    with open(expfile) as f:
        expected = f.readlines()
    assert mids.samfile_to_mids(infile) == expected
