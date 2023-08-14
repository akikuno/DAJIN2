from src.DAJIN2.core import report

###########################################################
# revcomp_cssplits
###########################################################


def test_revcomp_cssplits_substitution():
    cssplits = ["=A", "*AG", "=C", "=G"]
    test = report.report_mutation.revcomp_cssplits(cssplits)
    answer = ["=C", "=G", "*TC", "=T"]
    assert test == answer


def test_revcomp_cssplits_substitution_with_N():
    cssplits = ["=A", "*AN", "N", "=G"]
    test = report.report_mutation.revcomp_cssplits(cssplits)
    answer = ["=C", "N", "*TN", "=T"]
    assert test == answer


def test_revcomp_cssplits_insertion():
    cssplits = ["=A", "+A|+A|+T|+G|=A", "=G"]
    test = report.report_mutation.revcomp_cssplits(cssplits)
    answer = ["=C", "=T", "+C|+A|+T|+T|=T"]
    assert test == answer


def test_revcomp_cssplits_insertion_long():
    cssplits = ["=A", "=A", "+A|+A|+T|+G|=A", "=G", "=G"]
    test = report.report_mutation.revcomp_cssplits(cssplits)
    answer = ["=C", "=C", "=T", "+C|+A|+T|+T|=T", "=T"]
    assert test == answer


def test_revcomp_cssplits_insertion_with_substitution():
    cssplits = ["=A", "+A|+A|+T|+G|*AG", "=G"]
    test = report.report_mutation.revcomp_cssplits(cssplits)
    answer = ["=C", "*TC", "+C|+A|+T|+T|=T"]
    assert test == answer


def test_revcomp_cssplits_insertion_with_N():
    cssplits = ["N", "+A|+A|+T|+G|N", "=G"]
    test = report.report_mutation.revcomp_cssplits(cssplits)
    answer = ["=C", "N", "+C|+A|+T|+T|N"]
    assert test == answer


def test_revcomp_cssplits_deletion():
    cssplits = ["=A", "-C", "-A", "=G"]
    test = report.report_mutation.revcomp_cssplits(cssplits)
    answer = ["=C", "-T", "-G", "=T"]
    assert test == answer


def test_revcomp_cssplits_N():
    cssplits = ["=A", "N", "N", "=G"]
    test = report.report_mutation.revcomp_cssplits(cssplits)
    answer = ["=C", "N", "N", "=T"]
    assert test == answer


def test_revcomp_cssplits_inversion():
    cssplits = ["=A", "=a", "=g", "=G"]
    test = report.report_mutation.revcomp_cssplits(cssplits)
    answer = ["=C", "=c", "=t", "=T"]
    assert test == answer


###########################################################
# annotate inversion
###########################################################


def test_annotate_inversion():
    cssplits = ["=A", "=a", "=g", "=G"]
    test = report.report_mutation.annotate_inversion(cssplits)
    answer = ["=A", "@=a", "@=g", "=G"]
    assert test == answer


###########################################################
# group by mutation
###########################################################


def test_group_by_mutation_deletion():
    cssplits = ["=A", "-G", "-T", "=G"]
    test = report.report_mutation.group_by_mutation(cssplits)
    answer = [["=A"], ["-G", "-T"], ["=G"]]
    assert test == answer


def test_group_by_mutation_insertion():
    cssplits = ["=A", "+A|+A|=G", "=G"]
    test = report.report_mutation.group_by_mutation(cssplits)
    answer = [["=A"], ["+A|+A|=G"], ["=G"]]
    assert test == answer


def test_group_by_mutation_inversion():
    cssplits = ["=A", "@=a", "@=g", "=G"]
    test = report.report_mutation.group_by_mutation(cssplits)
    answer = [["=A"], ["@=a", "@=g"], ["=G"]]
    assert test == answer


########################################################################
# _report_mutations
########################################################################


def test_report_mutations_substitution():
    cssplits = ["=A", "*AG", "=G"]
    GENOME_COODINATES = {"genome": "mm10", "chrom": "chr1", "start": 0, "end": 2, "strand": "+"}
    header = "test"
    cssplits_inversion = report.report_mutation.annotate_inversion(cssplits)
    cssplits_grouped = report.report_mutation.group_by_mutation(cssplits_inversion)
    test = report.report_mutation.report_mutations(cssplits_grouped, GENOME_COODINATES, header)
    answer = [["test", "mm10", "chr1", 1, 1, "substitution: A>G"]]
    assert test == answer


def test_report_mutations_consecutive_substitution():
    cssplits = ["=A", "*AG", "*CT", "=G"]
    GENOME_COODINATES = {"genome": "mm10", "chrom": "chr1", "start": 0, "end": 3, "strand": "+"}
    header = "test"
    cssplits_inversion = report.report_mutation.annotate_inversion(cssplits)
    cssplits_grouped = report.report_mutation.group_by_mutation(cssplits_inversion)
    test = report.report_mutation.report_mutations(cssplits_grouped, GENOME_COODINATES, header)
    answer = [["test", "mm10", "chr1", 1, 1, "substitution: A>G"], ["test", "mm10", "chr1", 2, 2, "substitution: C>T"]]
    assert test == answer


def test_report_mutations_deletion():
    cssplits = ["=A", "-G", "-T", "=G"]
    GENOME_COODINATES = {"genome": "mm10", "chrom": "chr1", "start": 0, "end": 3, "strand": "+"}
    header = "test"
    cssplits_inversion = report.report_mutation.annotate_inversion(cssplits)
    cssplits_grouped = report.report_mutation.group_by_mutation(cssplits_inversion)
    test = report.report_mutation.report_mutations(cssplits_grouped, GENOME_COODINATES, header)
    answer = [["test", "mm10", "chr1", 1, 2, "2bp deletion: GT"]]
    assert test == answer


def test_report_mutations_insertion():
    cssplits = ["=A", "+G|+T|=A", "=G"]
    GENOME_COODINATES = {"genome": "mm10", "chrom": "chr1", "start": 0, "end": 3, "strand": "+"}
    header = "test"
    cssplits_inversion = report.report_mutation.annotate_inversion(cssplits)
    cssplits_grouped = report.report_mutation.group_by_mutation(cssplits_inversion)
    test = report.report_mutation.report_mutations(cssplits_grouped, GENOME_COODINATES, header)
    answer = [["test", "mm10", "chr1", 1, 1, "2bp insertion: GT"]]
    assert test == answer


def test_report_mutations_insertion_with_substitution():
    cssplits = ["=A", "+G|+T|*GA", "=G"]
    GENOME_COODINATES = {"genome": "mm10", "chrom": "chr1", "start": 0, "end": 3, "strand": "+"}
    header = "test"
    cssplits_inversion = report.report_mutation.annotate_inversion(cssplits)
    cssplits_grouped = report.report_mutation.group_by_mutation(cssplits_inversion)
    test = report.report_mutation.report_mutations(cssplits_grouped, GENOME_COODINATES, header)
    answer = [["test", "mm10", "chr1", 1, 1, "2bp insertion: GT"]]
    assert test == answer


def test_report_mutations_inversion():
    cssplits = ["=A", "=a", "=t", "=G"]
    GENOME_COODINATES = {"genome": "mm10", "chrom": "chr1", "start": 0, "end": 4, "strand": "+"}
    header = "test"
    cssplits_inversion = report.report_mutation.annotate_inversion(cssplits)
    cssplits_grouped = report.report_mutation.group_by_mutation(cssplits_inversion)
    test = report.report_mutation.report_mutations(cssplits_grouped, GENOME_COODINATES, header)
    answer = [["test", "mm10", "chr1", 1, 2, "2bp inversion: AT"]]
    assert test == answer


def test_report_mutations_various():
    cssplits = ["=A", "*AG", "+A|+A|-A", "N", "-G", "*CG", "*AG", "=a", "=t", "=A"]
    GENOME_COODINATES = {"genome": "mm10", "chrom": "chr1", "start": 0, "end": 4, "strand": "+"}
    header = "test"
    cssplits_inversion = report.report_mutation.annotate_inversion(cssplits)
    cssplits_grouped = report.report_mutation.group_by_mutation(cssplits_inversion)
    test = report.report_mutation.report_mutations(cssplits_grouped, GENOME_COODINATES, header)
    answer = [
        ["test", "mm10", "chr1", 1, 1, "substitution: A>G"],
        ["test", "mm10", "chr1", 2, 2, "2bp insertion: AA"],
        ["test", "mm10", "chr1", 2, 2, "1bp deletion: A"],
        ["test", "mm10", "chr1", 3, 3, "1bp unknown bases"],
        ["test", "mm10", "chr1", 4, 4, "1bp deletion: G"],
        ["test", "mm10", "chr1", 5, 5, "substitution: C>G"],
        ["test", "mm10", "chr1", 6, 6, "substitution: A>G"],
        ["test", "mm10", "chr1", 7, 8, "2bp inversion: AT"],
    ]
    assert test == answer


def test_report_mutations_genome_coodinates():
    cssplits = ["=A", "*AG", "=G"]
    GENOME_COODINATES = {"genome": "mm10", "chrom": "chrX", "start": 100, "end": 103, "strand": "+"}
    header = "test"
    cssplits_inversion = report.report_mutation.annotate_inversion(cssplits)
    cssplits_grouped = report.report_mutation.group_by_mutation(cssplits_inversion)
    test = report.report_mutation.report_mutations(cssplits_grouped, GENOME_COODINATES, header)
    answer = [["test", "mm10", "chrX", 101, 101, "substitution: A>G"]]
    assert test == answer
