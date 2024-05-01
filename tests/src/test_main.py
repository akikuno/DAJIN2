from __future__ import annotations

import pytest

from src.DAJIN2 import main


###########################################################
# validate_columns
###########################################################
def test_validate_columns_all_required_present():
    columns = {"sample", "control", "allele", "name"}
    # No error should be raised
    main.validate_headers_of_batch_file(columns, "test_filepath")


def test_validate_columns_missing_required():
    columns = {"sample", "control", "allele"}
    with pytest.raises(
        ValueError, match='test_filepath must contain "sample", "control", "allele" and "name" in the header'
    ):
        main.validate_headers_of_batch_file(columns, "test_filepath")


def test_validate_columns_extra_not_accepted():
    columns = {"sample", "control", "allele", "name", "extra_column"}
    with pytest.raises(
        ValueError,
        match='Accepted header names of test_filepath are "sample", "control", "allele", "name", or "genome".',
    ):
        main.validate_headers_of_batch_file(columns, "test_filepath")
