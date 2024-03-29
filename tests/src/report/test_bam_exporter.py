from __future__ import annotations

import pytest
from pathlib import Path
from unittest.mock import patch
from src.DAJIN2.core.report import bam_exporter

###############################################################################
# Output
###############################################################################


@pytest.fixture
def mock_pysam_sort_and_index():
    with patch("pysam.sort") as mock_sort:
        with patch("pysam.index") as mock_index:
            yield mock_sort, mock_index


def test_write_sam_to_bam(mock_pysam_sort_and_index):
    mock_sort, mock_index = mock_pysam_sort_and_index
    sam_data = [["r1", "seq1", "+", "chr1", "1000"], ["r2", "seq2", "-", "chr2", "2000"]]
    path_sam = "test.sam"
    path_bam = "test.bam"
    bam_exporter.write_sam_to_bam(sam_data, path_sam, path_bam)

    mock_sort.assert_called_once()
    mock_index.assert_called_once()

    with open(path_sam) as f:
        content = f.read()
        expected_content = "r1\tseq1\t+\tchr1\t1000\nr2\tseq2\t-\tchr2\t2000\n"
        assert content == expected_content

    # テスト後のクリーンアップ
    Path(path_sam).unlink()
