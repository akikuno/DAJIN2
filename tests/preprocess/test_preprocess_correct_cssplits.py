# """
# - repetitive dels
# - Shared mutation percentages in 5-mer
#     - sample:  [{"=A,=A,=A,=A,=A": [0.98,0.98,0.98,0.98,0.98,0.98]}, ]
#     - control: [{"=A,=A,=A,=A,=A": [0.98,0.98,0.98,0.98,0.98,0.98]}]
# """

# from misc.scripts import scratch_correct_cssplits

# ###############################################################################
# # extract_indexes_with_both_ends_not_N
# ###############################################################################


# def test_extract_indexes_with_both_ends_not_N___noN():
#     cssplits = ["=A,=C,=G,=T"]
#     cssplits = [cs.split(",") for cs in cssplits]
#     test = scratch_correct_cssplits.extract_indexes_with_both_ends_not_N(cssplits)
#     answer = [(0, 3)]
#     assert test == answer


# def test_extract_indexes_with_both_ends_not_N___leftN():
#     cssplits = ["N,N,N,N,=A,=C,=G,=T"]
#     cssplits = [cs.split(",") for cs in cssplits]
#     test = scratch_correct_cssplits.extract_indexes_with_both_ends_not_N(cssplits)
#     answer = [(4, 7)]
#     assert test == answer


# def test_extract_indexes_with_both_ends_not_N___rightN():
#     cssplits = ["=A,=C,=G,=T,N,N,N,N"]
#     cssplits = [cs.split(",") for cs in cssplits]
#     test = scratch_correct_cssplits.extract_indexes_with_both_ends_not_N(cssplits)
#     answer = [(0, 3)]
#     assert test == answer


# def test_extract_indexes_with_both_ends_not_N___bothN():
#     cssplits = ["N,N,=A,=C,=G,=T,=A,N,N"]
#     cssplits = [cs.split(",") for cs in cssplits]
#     test = scratch_correct_cssplits.extract_indexes_with_both_ends_not_N(cssplits)
#     answer = [(2, 6)]
#     assert test == answer


# def test_extract_indexes_with_both_ends_not_N___multiple():
#     cssplits = ["=A,=C,=G,=T", "N,N,N,N,=A,=C,=G,=T", "=A,=C,=G,=T,N,N,N,N", "N,N,=A,=C,=G,=T,=A,N,N"]
#     cssplits = [cs.split(",") for cs in cssplits]
#     test = scratch_correct_cssplits.extract_indexes_with_both_ends_not_N(cssplits)
#     answer = [(0, 3), (4, 7), (0, 3), (2, 6)]
#     assert test == answer


# ###############################################################################
# # call_count
# ###############################################################################


# def test_call_count():
#     cssplits = ["=A,+G|=C,-G,*TA,=A", "=A,+G|=C,-G,=T,=A"]
#     cssplits = [cs.split(",") for cs in cssplits]
#     indexes = [(0, 4), (0, 4)]
#     test = scratch_correct_cssplits.call_count(cssplits, indexes)
#     answer = {1: {"=A,+G|=C,-G": 2}, 2: {"+G|=C,-G,*TA": 1, "+G|=C,-G,=T": 1}, 3: {"-G,*TA,=A": 1}}
#     assert test == answer


# def test_call_count_N():
#     cssplits = ["=A,+G|=C,-G,*TA,=A", "N,N,-G,*TA,=A"]
#     cssplits = [cs.split(",") for cs in cssplits]
#     indexes = [(0, 4), (2, 4)]
#     test = scratch_correct_cssplits.call_count(cssplits, indexes)
#     answer = {1: {"=A,+G|=C,-G": 1}, 2: {"+G|=C,-G,*TA": 1}, 3: {"-G,*TA,=A": 2}}
#     assert test == answer


# ###############################################################################
# # call_percentage
# ###############################################################################


# def test_call_percentage():
#     cssplits = ["=A,+G|=C,-G,*TA,=A", "=A,+G|=C,-G,=T,=A"]
#     cssplits = [cs.split(",") for cs in cssplits]
#     counts = {1: {"=A,+G|=C,-G": 2}, 2: {"+G|=C,-G,*TA": 1, "+G|=C,-G,=T": 1}, 3: {"-G,*TA,=A": 1}}
#     test = scratch_correct_cssplits.call_percentage(cssplits, counts)
#     answer = {1: {"=A,+G|=C,-G": 100.0}, 2: {"+G|=C,-G,*TA": 50.0, "+G|=C,-G,=T": 50.0}, 3: {"-G,*TA,=A": 50.0}}
#     assert test == answer


# ###############################################################################
# # subtract_percentage
# ###############################################################################


# def test_subtract_percentage():
#     percent_sample = {
#         1: {"=A,+G|=C,-G": 100.0},
#         2: {"+G|=C,-G,*TA": 50.0, "+G|=C,-G,=T": 50.0},
#         3: {"-G,*TA,=A": 50.0, "-G,-T,=A": 50.0},
#     }
#     percent_control = {
#         1: {"=A,+G|=C,-G": 75.0},
#         2: {"+G|=C,-G,*TA": 50.0, "+G|=C,-G,=T": 50.0},
#         3: {"-G,*TA,=A": 25.0},
#     }
#     test = scratch_correct_cssplits.subtract_percentage(percent_sample, percent_control)
#     answer = {
#         1: {"=A,+G|=C,-G": 25.0},
#         2: {"+G|=C,-G,*TA": 0.0, "+G|=C,-G,=T": 0.0},
#         3: {"-G,*TA,=A": 25.0, "-G,-T,=A": 50.0},
#     }
#     assert test == answer


# ###############################################################################
# # select_candidate_mutation
# ###############################################################################


# def test_select_candidate_mutation():
#     percent_subtracted = {
#         1: {"=A,+G|=C,-G": 25.0},
#         2: {"+G|=C,-G,*TA": 0.0, "+G|=C,-G,=T": 0.0},
#         3: {"-G,*TA,=A": 25.0, "-G,-T,=A": 50.0},
#     }
#     test = scratch_correct_cssplits.select_candidate_mutation(percent_subtracted)
#     answer = {1: {"=A,+G|=C,-G"}, 2: set(), 3: {"-G,*TA,=A", "-G,-T,=A"}}
#     assert test == answer


# ###############################################################################
# # update_cssplits
# ###############################################################################


# def test_update_cssplits():
#     cssplits = ["=A,+G|=C,-G,*TA,=A", "=A,+G|=C,-G,-T,=A"]
#     cssplits = [cs.split(",") for cs in cssplits]
#     sequence = "ACGTA"
#     candidate_mutation = {1: {"=A,+G|=C,-G"}, 2: set(), 3: {"-G,*TA,=A", "-G,-T,=A"}}
#     test = scratch_correct_cssplits.update_cssplits(cssplits, sequence, candidate_mutation)
#     answer = [["=A", "+G|=C", "=G", "*TA", "=A"], ["=A", "+G|=C", "=G", "-T", "=A"]]
#     assert test == answer

