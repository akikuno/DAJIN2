# tests/test_bed_utils.py
# Pytest unit tests for parse_bed_file and bed_to_genome_coordinates.
# Assumes the module under test is saved as "bed_utils.py" in the test root.
# If your filename is different, change the import line below accordingly.

from pathlib import Path

import pytest

from DAJIN2.utils.bed_handler import BEDError, bed_to_genome_coordinates, parse_bed_file


def write_text(tmp_path: Path, name: str, content: str) -> Path:
    p = tmp_path / name
    p.write_text(content, encoding="utf-8")
    return p


@pytest.mark.parametrize(
    "bed_text,expected",
    [
        (
            # Valid single BED6 line, plus comments/track lines that should be ignored
            "\n".join(
                [
                    "# comment",
                    "track name=mybed",
                    "chr1\t0\t100\tmm39\t248956422\t+",
                ]
            ),
            [
                {
                    "chrom": "chr1",
                    "start": 0,
                    "end": 100,
                    "genome": "mm39",
                    "chrom_size": 248956422,
                    "strand": "+",
                }
            ],
        ),
        (
            # Valid minus-strand
            "chrX\t10\t20\tmm39\t156040895\t-",
            [
                {
                    "chrom": "chrX",
                    "start": 10,
                    "end": 20,
                    "genome": "mm39",
                    "chrom_size": 156040895,
                    "strand": "-",
                }
            ],
        ),
    ],
)
def test_parse_bed_file_valid_cases(tmp_path, bed_text, expected):
    bed = write_text(tmp_path, "ok.bed", bed_text)
    intervals = parse_bed_file(bed)
    assert intervals == expected


@pytest.mark.parametrize(
    "bed_text,err_substr",
    [
        # Wrong number of fields
        ("chr1\t0\t100", "Invalid BED format"),
        ("chr1\t0\t100\tmm39\t248956422", "Invalid BED format"),
        ("chr1\t0\t100\tmm39\t248956422\t+\textra", "Invalid BED format"),
        # Missing/invalid strand
        ("chr1\t0\t100\tmm39\t248956422\t*", "Invalid or missing strand"),
        ("chr1\t0\t100\tmm39\t248956422\t.", "Invalid or missing strand"),
        # Coordinate validation
        ("chr1\t-1\t10\tmm39\t248956422\t+", "Invalid start position"),
        ("chr1\t10\t10\tmm39\t248956422\t+", "Invalid end position"),
        ("chr1\tA\t10\tmm39\t248956422\t+", "Invalid coordinate format"),
        ("chr1\t10\tB\tmm39\t248956422\t+", "Invalid coordinate format"),
        # Chrom size validation
        ("chr1\t0\t10\tmm39\tNA\t+", "Invalid chromosome size format"),
        ("chr1\t0\t10\tmm39\t0\t+", "Invalid chromosome size"),
        ("chr1\t0\t10\tmm39\t-5\t+", "Invalid chromosome size"),
    ],
)
def test_parse_bed_file_invalid_cases(tmp_path, bed_text, err_substr):
    bed = write_text(tmp_path, "bad.bed", bed_text)
    with pytest.raises(BEDError) as ei:
        parse_bed_file(bed)
    assert err_substr in str(ei.value)


def test_parse_bed_file_empty_or_comments_only(tmp_path):
    bed = write_text(tmp_path, "empty.bed", "# only comments\ntrack name=x\n")
    with pytest.raises(BEDError) as ei:
        parse_bed_file(bed)
    assert "No valid intervals found" in str(ei.value)


def test_parse_bed_file_file_not_found(tmp_path):
    missing = tmp_path / "nope.bed"
    with pytest.raises(FileNotFoundError) as ei:
        parse_bed_file(missing)
    assert "BED file not found" in str(ei.value)


def test_bed_to_genome_coordinates_single_interval(tmp_path):
    bed = write_text(tmp_path, "one.bed", "chr2\t100\t200\thg38\t242193529\t+")
    got = bed_to_genome_coordinates(bed)
    assert got == {
        "genome": "hg38",
        "chrom": "chr2",
        "start": 100,
        "end": 200,
        "strand": "+",
        "chrom_size": 242193529,
    }


def test_bed_to_genome_coordinates_multiple_intervals_warn_and_use_first(tmp_path, caplog):
    bed = write_text(
        tmp_path,
        "multi.bed",
        "\n".join(
            [
                "chr3\t10\t20\tmm39\t156040895\t+",
                "chr4\t30\t40\tmm39\t156040895\t-",
            ]
        ),
    )
    with caplog.at_level("WARNING"):
        got = bed_to_genome_coordinates(bed)
        # Check warning message mentions count and first interval
        assert any(
            "BED file contains 2 intervals. Using first interval: chr3:10-20" in r.message for r in caplog.records
        )
    assert got["chrom"] == "chr3"
    assert got["start"] == 10
    assert got["end"] == 20
    assert got["strand"] == "+"
    assert got["genome"] == "mm39"
    assert got["chrom_size"] == 156040895


def test_parse_bed_file_integrity_double_validation(tmp_path):
    # start < end validation is checked twice; ensure it trips properly
    bed = write_text(tmp_path, "bad_coords.bed", "chr5\t50\t50\tmm39\t156040895\t+")
    with pytest.raises(BEDError) as ei:
        parse_bed_file(bed)
    msg = str(ei.value)
    assert "Invalid end position" in msg or "must be > start" in msg
