import midsv
from pathlib import Path
from importlib import reload
from src.DAJIN2.report import report_bam

reload(report_bam)


def test_remove_microhomology():
    sam = midsv.read_sam("tests/data/report_bam/microhomology-deletion.sam")
    test = report_bam.remove_microhomology(sam)
    answer = Path("tests/data/report_bam/answer.txt").read_text()
    answer = eval(answer)
    assert test == answer


def test_remove_microhomology_real_singe_read():
    seq_id = "0ef147016ef0"
    sam = midsv.read_sam("tests/data/report_bam/barcode54_allele2_before.sam")
    sam_headers = [s for s in sam if s[0].startswith("@")]
    sam_contents = [s for s in sam if seq_id in s[0] and s[9] != "*"]
    sam = sam_headers + sam_contents
    test = report_bam.remove_microhomology(sam)
    answer = midsv.read_sam(f"tests/data/report_bam/barcode54_allele2_{seq_id}_after.sam")
    assert test == answer


def test_remove_microhomology_real_all_500_reads():
    sam = midsv.read_sam("tests/data/report_bam/barcode54_allele2_before.sam")
    test = report_bam.remove_microhomology(sam)
    answer = midsv.read_sam("tests/data/report_bam/barcode54_allele2_after.sam")
    assert test == answer

