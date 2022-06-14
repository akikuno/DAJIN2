from src.DAJIN2.utils.io import *
from pathlib import Path
import tempfile
import filecmp

input_path = str(Path("tests", "io", "exam.jsonl"))
output_path = tempfile.NamedTemporaryFile().name

exam = [{"name": "nobita", "math": 2, "english": 5}, {"name": "shizuka", "math": 100, "english": 95}]


def test_read_jsonl():
    exam_eval = read_jsonl(input_path)
    assert exam == exam_eval


def test_write_jsonl():
    write_jsonl(exam, output_path)
    assert filecmp.cmp(input_path, output_path)

