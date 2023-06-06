import midsv
from pathlib import Path
from importlib import reload
from src.DAJIN2.core.report import report_bam
from src.DAJIN2.core.report.remove_microhomology import remove_microhomology

reload(report_bam)


def test_remove_overlapped_reads():
    """
    # Overlapped reads in .tmpDAJIN/sam/barcode31_control.sam:
        - 01f0e046-b6dc-46f3-ad5b-96ed91484778
    # Non-overlapped reads in .tmpDAJIN/sam/barcode31_control.sam
        - a224e9ca-d634-4490-bf77-abd9d5cbd7bc
    """
    sam = midsv.read_sam("tests/data/report/report_bam/remove_overlap.sam")
    sam = list(sam)
    test = report_bam.remove_overlapped_reads(sam)
    answer = midsv.read_sam("tests/data/report/report_bam/answer.sam")
    answer = list(answer)
    assert test == answer


def test_remove_microhomology():
    sam = midsv.read_sam("tests/data/report/report_bam/microhomology-deletion.sam")
    sam = list(sam)
    test = remove_microhomology(sam)
    answer = Path("tests/data/report/report_bam/answer.txt").read_text()
    answer = eval(answer)
    assert test == answer


def test_remove_microhomology_overlapped_softclip():
    sam = midsv.read_sam("tests/data/report/report_bam/cables2_barcode34_microhomology.sam")
    sam = list(sam)
    test = remove_microhomology(sam)
    answer = list(midsv.read_sam("tests/data/report/report_bam/answer_cables2_barcode34_microhomology.sam"))
    assert test == answer


def test_remove_microhomology_real_singe_read():
    seq_id = "0ef147016ef0"
    sam = midsv.read_sam("tests/data/report/report_bam/barcode54_allele2_before.sam")
    sam = list(sam)
    sam_headers = [s for s in sam if s[0].startswith("@")]
    sam_contents = [s for s in sam if seq_id in s[0] and s[9] != "*"]
    sam = sam_headers + sam_contents
    test = remove_microhomology(sam)
    answer = midsv.read_sam(f"tests/data/report/report_bam/barcode54_allele2_{seq_id}_after.sam")
    answer = list(answer)
    assert test == answer


def test_remove_microhomology_real_all_500_reads():
    sam = midsv.read_sam("tests/data/report/report_bam/barcode54_allele2_before.sam")
    sam = list(sam)
    test = remove_microhomology(sam)
    answer = midsv.read_sam("tests/data/report/report_bam/barcode54_allele2_after.sam")
    answer = list(answer)
    assert test == answer
