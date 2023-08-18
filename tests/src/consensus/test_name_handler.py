from collections import defaultdict
from DAJIN2.core.consensus.name_handler import (
    _detect_sv,
    call_allele_name,
    update_key_by_allele_name,
    add_key_by_allele_name,
)

###########################################################
# detect_sv
###########################################################


def test_detect_sv_threshold():
    cons_percentages = defaultdict(list)
    cons_percentages["sample1"] = [{"N": 90, "=C": 10}]  # one n
    cons_percentages["sample2"] = [{"+G|+G|=A": 80, "=C": 20}]  # two insertion
    cons_percentages["sample3"] = [{"-A": 100}, {"-A": 100}, {"-A": 100}]  # three deletion
    cons_percentages["sample4"] = [{"*AT": 100}, {"*AT": 100}, {"*AT": 100}, {"*AT": 100}]  # four substitution
    cons_percentages["sample5"] = [{"=a": 100}]  # inversion
    assert _detect_sv(cons_percentages, threshold=1) == [True, True, True, True, True]
    assert _detect_sv(cons_percentages, threshold=2) == [False, True, True, True, True]
    assert _detect_sv(cons_percentages, threshold=3) == [False, False, True, True, True]
    assert _detect_sv(cons_percentages, threshold=4) == [False, False, False, True, True]
    assert _detect_sv(cons_percentages, threshold=5) == [False, False, False, False, True]
    assert _detect_sv(cons_percentages, threshold=6) == [False, False, False, False, True]


###########################################################
# call_allele_name
###########################################################


def test_call_allele_name():
    pass


#     cons_sequences = defaultdict(dict)
#     cons_sequences[(1, 10, 90)] = "ATGC"
#     cons_sequences[(2, 20, 80)] = "CGTA"

#     cons_percentages = defaultdict(list)
#     cons_percentages[10] = [{"A": 30, "T": 20}, {"N": 40, "G": 60}]
#     cons_percentages[20] = [{"A": 40, "T": 60}, {"+G|": 55, "C": 45}]

#     FASTA_ALLELES = {1: "ATGC", 2: "GTCA"}

#     expected = {10: "allele10_1_sv_90%", 20: "allele20_2_sv_80%"}
#     assert call_allele_name(cons_sequences, cons_percentages, FASTA_ALLELES) == expected


def test_update_key_by_allele_name():
    pass


#     cons = {10: "sample_data_1", 20: "sample_data_2"}
#     allele_names = {10: "allele10_1_sv_90%", 20: "allele20_2_sv_80%"}

#     expected = {"allele10_1_sv_90%": "sample_data_1", "allele20_2_sv_80%": "sample_data_2"}
#     assert update_key_by_allele_name(cons, allele_names) == expected


def test_add_key_by_allele_name():
    pass


#     clust_sample = [{"LABEL": 10, "DATA": "sample1"}, {"LABEL": 20, "DATA": "sample2"}]
#     allele_names = {10: "allele10_1_sv_90%", 20: "allele20_2_sv_80%"}

#     expected = [
#         {"LABEL": 10, "DATA": "sample1", "NAME": "allele10_1_sv_90%"},
#         {"LABEL": 20, "DATA": "sample2", "NAME": "allele20_2_sv_80%"},
#     ]
#     assert add_key_by_allele_name(clust_sample, allele_names) == expected
