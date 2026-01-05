from __future__ import annotations

import gzip
from pathlib import Path
from tempfile import TemporaryDirectory

import pytest

from DAJIN2.utils import fastx_handler

########################################################################
# Extract basename
########################################################################


def test_extract_filename():
    fastq_path = "tests/hoge/fuga/example.fq.gz"
    test = fastx_handler.extract_filename(fastq_path)
    answer = "example"
    assert test == answer


def test_extract_filename_fail():
    fastq_path = "tests/hoge/fuga/.fq.gz"
    with pytest.raises(ValueError) as e:
        fastx_handler.extract_filename(fastq_path)
    assert str(e.value) == "Provided name is empty or consists only of whitespace"


def test_extract_filename_change_filename():
    fastq_path = "tests/hoge/fuga/x?y|z.fq.gz"
    test = fastx_handler.extract_filename(fastq_path)
    answer = "x-y-z"
    assert test == answer


########################################################################
# Convert allele file to dictionary type fasta format
########################################################################


def test_dictionize_allele_fasta():
    path_fasta = "tests/data/preprocess/format_input/fasta.fa"
    test = fastx_handler.dictionize_allele(path_fasta)
    answer = {"test1": "ACGTACGTACGTACGTACGTACGTACGTACGTACGT", "test2": "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"}
    assert test == answer


def test_dictionize_allele_fasta_wrap():
    path_fasta = "tests/data/preprocess/format_input/fasta_wrap.fa"
    test = fastx_handler.dictionize_allele(path_fasta)
    answer = {"test1": "ACGTACGTACGTACGTACGTACGTACGTACGTACGT", "test2": "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"}
    assert test == answer


########################################################################
# convert_fasta_to_fastq
########################################################################


@pytest.mark.parametrize(
    "fasta_content, expected_fastq_content",
    [
        (
            ">seq1\nAGCTTAGCTAG\n>seq2\nCGATCGTAGC\n",
            "@seq1\nAGCTTAGCTAG\n+\nIIIIIIIIIII\n@seq2\nCGATCGTAGC\n+\nIIIIIIIIII\n",
        ),
    ],
)
def test_convert_fasta_to_fastq(tmp_path, fasta_content, expected_fastq_content):
    fasta_path = tmp_path / "test.fasta"
    fastq_path = tmp_path / "test.fastq.gz"

    # Write test FASTA content
    fasta_path.write_text(fasta_content)

    # Run conversion function
    fastx_handler.convert_fasta_to_fastq(fasta_path, fastq_path)

    # Read and check if gzipped FASTQ output matches expected output
    with gzip.open(fastq_path, "rt") as f:
        result = f.read()
    assert result == expected_fastq_content


@pytest.mark.parametrize(
    "filenames, suffixes, expected",
    [
        # 単純なケース：suffix一致
        (["a.txt", "b.csv", "c.txt"], {".txt"}, ["a.txt", "c.txt"]),
        # .fai, .baiを除外
        (["ref.fai", "sample.bam", "sample.bai"], {".bam"}, ["sample.bam"]),
        # 一致なし
        (["a.txt", "b.csv"], {".tsv"}, []),
        # 複数のsuffixを許可
        (["a.txt", "b.csv", "c.tsv"], {".csv", ".tsv"}, ["b.csv", "c.tsv"]),
        # ディレクトリが空の場合
        ([], {".txt"}, []),
    ],
)
def test_get_path_input_files(filenames, suffixes, expected):
    with TemporaryDirectory() as tmpdir:
        dir_path = Path(tmpdir)
        for name in filenames:
            (dir_path / name).touch()
        result = fastx_handler._get_path_input_files(dir_path, suffixes)
        result_names = sorted([p.name for p in result])
        assert result_names == sorted(expected)
