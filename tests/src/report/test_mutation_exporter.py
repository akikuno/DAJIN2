from __future__ import annotations

from src.DAJIN2.core import report

###########################################################
# revcomp_midsvs
###########################################################


def test_revcomp_midsvs_substitution():
    midsv_tags = ["=A", "*AG", "=C", "=G"]
    test = report.mutation_exporter.revcomp_midsvs(midsv_tags)
    answer = ["=C", "=G", "*TC", "=T"]
    assert test == answer


def test_revcomp_midsvs_substitution_with_N():
    midsv_tags = ["=A", "*AN", "=N", "=G"]
    test = report.mutation_exporter.revcomp_midsvs(midsv_tags)
    answer = ["=C", "=N", "*TN", "=T"]
    assert test == answer


def test_revcomp_midsvs_insertion():
    midsv_tags = ["=A", "+A|+A|+T|+G|=A", "=G"]
    test = report.mutation_exporter.revcomp_midsvs(midsv_tags)
    answer = ["=C", "=T", "+C|+A|+T|+T|=T"]
    assert test == answer


def test_revcomp_midsvs_insertion_long():
    midsv_tags = ["=A", "=A", "+A|+A|+T|+G|=A", "=G", "=G"]
    test = report.mutation_exporter.revcomp_midsvs(midsv_tags)
    answer = ["=C", "=C", "=T", "+C|+A|+T|+T|=T", "=T"]
    assert test == answer


def test_revcomp_midsvs_insertion_with_substitution():
    midsv_tags = ["=A", "+A|+A|+T|+G|*AG", "=G"]
    test = report.mutation_exporter.revcomp_midsvs(midsv_tags)
    answer = ["=C", "*TC", "+C|+A|+T|+T|=T"]
    assert test == answer


def test_revcomp_midsvs_insertion_with_N():
    midsv_tags = ["=N", "+A|+A|+T|+G|=N", "=G"]
    test = report.mutation_exporter.revcomp_midsvs(midsv_tags)
    answer = ["=C", "=N", "+C|+A|+T|+T|=N"]
    assert test == answer


def test_revcomp_midsvs_deletion():
    midsv_tags = ["=A", "-C", "-A", "=G"]
    test = report.mutation_exporter.revcomp_midsvs(midsv_tags)
    answer = ["=C", "-T", "-G", "=T"]
    assert test == answer


def test_revcomp_midsvs_N():
    midsv_tags = ["=A", "=N", "=N", "=G"]
    test = report.mutation_exporter.revcomp_midsvs(midsv_tags)
    answer = ["=C", "=N", "=N", "=T"]
    assert test == answer


def test_revcomp_midsvs_inversion():
    midsv_tags = ["=A", "=a", "=g", "=G"]
    test = report.mutation_exporter.revcomp_midsvs(midsv_tags)
    answer = ["=C", "=c", "=t", "=T"]
    assert test == answer


###########################################################
# annotate inversion
###########################################################


def test_annotate_inversion():
    midsv_tags = ["=A", "=a", "=g", "=G"]
    test = report.mutation_exporter.annotate_inversion(midsv_tags)
    answer = ["=A", "@=a", "@=g", "=G"]
    assert test == answer


###########################################################
# group by mutation
###########################################################


def test_group_by_mutation_deletion():
    midsv_tags = ["=A", "-G", "-T", "=G"]
    test = report.mutation_exporter.group_by_mutation(midsv_tags)
    answer = [["=A"], ["-G", "-T"], ["=G"]]
    assert test == answer


def test_group_by_mutation_insertion():
    midsv_tags = ["=A", "+A|+A|=G", "=G"]
    test = report.mutation_exporter.group_by_mutation(midsv_tags)
    answer = [["=A"], ["+A|+A|=G"], ["=G"]]
    assert test == answer


def test_group_by_mutation_inversion():
    midsv_tags = ["=A", "@=a", "@=g", "=G"]
    test = report.mutation_exporter.group_by_mutation(midsv_tags)
    answer = [["=A"], ["@=a", "@=g"], ["=G"]]
    assert test == answer


########################################################################
# _report_mutations
########################################################################


def test_report_mutations_substitution():
    midsv_tags = ["=A", "*AG", "=G"]
    GENOME_COODINATES = {"genome": "mm10", "chrom": "chr1", "start": 0, "end": 2, "strand": "+"}
    header = "test"
    midsv_tags_inversion = report.mutation_exporter.annotate_inversion(midsv_tags)
    midsv_tags_grouped = report.mutation_exporter.group_by_mutation(midsv_tags_inversion)
    test = report.mutation_exporter.report_mutations(midsv_tags_grouped, GENOME_COODINATES, header)
    answer = [["test", "mm10", "chr1", "1", "1", "substitution: A>G"]]
    assert test == answer


def test_report_mutations_consecutive_substitution():
    midsv_tags = ["=A", "*AG", "*CT", "=G"]
    GENOME_COODINATES = {"genome": "mm10", "chrom": "chr1", "start": 0, "end": 3, "strand": "+"}
    header = "test"
    midsv_tags_inversion = report.mutation_exporter.annotate_inversion(midsv_tags)
    midsv_tags_grouped = report.mutation_exporter.group_by_mutation(midsv_tags_inversion)
    test = report.mutation_exporter.report_mutations(midsv_tags_grouped, GENOME_COODINATES, header)
    answer = [
        ["test", "mm10", "chr1", "1", "1", "substitution: A>G"],
        ["test", "mm10", "chr1", "2", "2", "substitution: C>T"],
    ]
    assert test == answer


def test_report_mutations_deletion():
    midsv_tags = ["=A", "-G", "-T", "=G"]
    GENOME_COODINATES = {"genome": "mm10", "chrom": "chr1", "start": 0, "end": 3, "strand": "+"}
    header = "test"
    midsv_tags_inversion = report.mutation_exporter.annotate_inversion(midsv_tags)
    midsv_tags_grouped = report.mutation_exporter.group_by_mutation(midsv_tags_inversion)
    test = report.mutation_exporter.report_mutations(midsv_tags_grouped, GENOME_COODINATES, header)
    answer = [["test", "mm10", "chr1", "1", "2", "2bp deletion: GT"]]
    assert test == answer


def test_report_mutations_insertion():
    midsv_tags = ["=A", "+G|+T|=A", "=G"]
    GENOME_COODINATES = {"genome": "mm10", "chrom": "chr1", "start": 0, "end": 3, "strand": "+"}
    header = "test"
    midsv_tags_inversion = report.mutation_exporter.annotate_inversion(midsv_tags)
    midsv_tags_grouped = report.mutation_exporter.group_by_mutation(midsv_tags_inversion)
    test = report.mutation_exporter.report_mutations(midsv_tags_grouped, GENOME_COODINATES, header)
    answer = [["test", "mm10", "chr1", "1", "1", "2bp insertion: GT"]]
    assert test == answer


def test_report_mutations_insertion_with_substitution():
    midsv_tags = ["=A", "+G|+T|*GA", "=G"]
    GENOME_COODINATES = {"genome": "mm10", "chrom": "chr1", "start": 0, "end": 3, "strand": "+"}
    header = "test"
    midsv_tags_inversion = report.mutation_exporter.annotate_inversion(midsv_tags)
    midsv_tags_grouped = report.mutation_exporter.group_by_mutation(midsv_tags_inversion)
    test = report.mutation_exporter.report_mutations(midsv_tags_grouped, GENOME_COODINATES, header)
    answer = [["test", "mm10", "chr1", "1", "1", "2bp insertion: GT"]]
    assert test == answer


def test_report_mutations_inversion():
    midsv_tags = ["=A", "=a", "=t", "=G"]
    GENOME_COODINATES = {"genome": "mm10", "chrom": "chr1", "start": 0, "end": 4, "strand": "+"}
    header = "test"
    midsv_tags_inversion = report.mutation_exporter.annotate_inversion(midsv_tags)
    midsv_tags_grouped = report.mutation_exporter.group_by_mutation(midsv_tags_inversion)
    test = report.mutation_exporter.report_mutations(midsv_tags_grouped, GENOME_COODINATES, header)
    answer = [["test", "mm10", "chr1", "1", "2", "2bp inversion: AT"]]
    assert test == answer


def test_report_mutations_various():
    midsv_tags = ["=A", "*AG", "+A|+A|-A", "=N", "-G", "*CG", "*AG", "=a", "=t", "=A"]
    GENOME_COODINATES = {"genome": "mm10", "chrom": "chr1", "start": 0, "end": 4, "strand": "+"}
    header = "test"
    midsv_tags_inversion = report.mutation_exporter.annotate_inversion(midsv_tags)
    midsv_tags_grouped = report.mutation_exporter.group_by_mutation(midsv_tags_inversion)
    test = report.mutation_exporter.report_mutations(midsv_tags_grouped, GENOME_COODINATES, header)
    answer = [
        ["test", "mm10", "chr1", "1", "1", "substitution: A>G"],
        ["test", "mm10", "chr1", "2", "2", "2bp insertion: AA"],
        ["test", "mm10", "chr1", "2", "2", "1bp deletion: A"],
        ["test", "mm10", "chr1", "4", "4", "1bp deletion: G"],
        ["test", "mm10", "chr1", "5", "5", "substitution: C>G"],
        ["test", "mm10", "chr1", "6", "6", "substitution: A>G"],
        ["test", "mm10", "chr1", "7", "8", "2bp inversion: AT"],
    ]
    assert test == answer


def test_report_mutations_genome_coodinates():
    midsv_tags = ["=A", "*AG", "=G"]
    GENOME_COODINATES = {"genome": "mm10", "chrom": "chrX", "start": 100, "end": 103, "strand": "+"}
    header = "test"
    midsv_tags_inversion = report.mutation_exporter.annotate_inversion(midsv_tags)
    midsv_tags_grouped = report.mutation_exporter.group_by_mutation(midsv_tags_inversion)
    test = report.mutation_exporter.report_mutations(midsv_tags_grouped, GENOME_COODINATES, header)
    answer = [["test", "mm10", "chrX", "101", "101", "substitution: A>G"]]
    assert test == answer
