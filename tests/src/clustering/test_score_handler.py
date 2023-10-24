from DAJIN2.core.clustering import score_handler


def test_call_count():
    cssplits_sample = [["=A", "*ag", "-gg", "=T"], ["=A", "-a", "-gg", "*ag"]]
    result = score_handler.call_count(iter(cssplits_sample))
    expected_result = [
        {"=A": 2},
        {"*ag": 1, "-a": 1},
        {"-gg": 2},
        {"=T": 1, "*ag": 1},
    ]
    assert result == expected_result


def test_call_percentage():
    counts = [
        {"=A": 2},
        {"*ag": 1, "-a": 1},
        {"-gg": 2},
        {"=T": 1, "*ag": 1},
    ]
    result = score_handler.call_percent(counts)
    expected_result = [
        {"=A": 100},
        {"*ag": 50, "-a": 50},
        {"-gg": 100},
        {"=T": 50, "*ag": 50},
    ]
    assert result == expected_result


def test_subtract_percentage():
    percent_sample = [{"A": 40, "T": 60}, {"C": 40, "G": 60}, {"A": 30, "T": 70}]
    percent_control = [{"A": 30, "T": 50}, {"C": 30, "G": 70}, {"A": 20, "T": 80}]
    knockin_loci = {1}
    result = score_handler.subtract_percentage(percent_sample, percent_control, knockin_loci)
    expected_result = [{"A": 10, "T": 10}, {"C": 40, "G": 60}, {"A": 10}]
    assert result == expected_result


def test_discard_common_error():
    percent_subtracted = [{"A": 0.4, "T": 99.6}, {"C": 0.5, "G": 99.5}, {"A": 0.6, "T": 99.4}]
    result = score_handler.discard_common_error(percent_subtracted, threshold=0.5)
    expected_result = [{"T": 99.6}, {"G": 99.5}, {"A": 0.6, "T": 99.4}]
    assert result == expected_result


def test_discard_matches_and_ns():
    percent_discarded = [
        {"=A,=A,=A": 10, "=A,+I,=A": 20},
        {"=C,=C,=C": 15, "=C,+I,=C": 25},
        {"=A,N,=A": 10, "=A,+I,=A": 10},
        {"=A,=A,=A": 10, "=A,N,=A": 40},
        {},
    ]
    result = score_handler.discard_matches_and_ns(percent_discarded)
    expected_result = [{"=A,+I,=A": 20}, {"=C,+I,=C": 25}, {"=A,+I,=A": 10}, {}, {}]
    assert result == expected_result


###############################################################################
# Handling insertions
###############################################################################


def test_group_consecutive_insertions():
    mutation_loci = [set(), set("+"), set("*"), set("+"), set("+"), set("+"), set(), set("+")]
    result = score_handler.group_consecutive_insertions(mutation_loci)
    assert result == [(1,), (3, 4, 5), (7,)]


def test_find_max_insertion_score_in_group():
    percent_discarded = [
        {"=A,=A,=A": 10, "=A,+I,=A": 20},
        {"=C,=C,=C": 15, "=C,+I,=C": 25},
        {"=A,=A,=A": 10, "=A,+I,=A": 10},
        {"=A,=A,=A": 10, "=A,+I,=A": 40},
    ]
    index_group = (2, 3)
    result = score_handler.find_max_insertion_score_in_group(percent_discarded, index_group)
    expected_result = 40
    assert result == expected_result


def test_update_insertion_scores_in_group_to_max():
    percent_discarded = [
        {"=A,=A,=A": 10, "=A,+I,=A": 20},
        {"=C,=C,=C": 15, "=C,+I,=C": 25},
        {"=A,=A,=A": 10, "=A,+I,=A": 10},
        {"=A,=A,=A": 10, "=A,+I,=A": 40},
    ]
    index_group = (2, 3)
    max_score = 40
    result = score_handler.update_insertion_scores_in_group_to_max(percent_discarded, index_group, max_score)
    expected_result = [
        {"=A,=A,=A": 10, "=A,+I,=A": 20},
        {"=C,=C,=C": 15, "=C,+I,=C": 25},
        {"=A,=A,=A": 10, "=A,+I,=A": 40},
        {"=A,=A,=A": 10, "=A,+I,=A": 40},
    ]
    assert result == expected_result


def test_update_insertion_score():
    percent_discarded = [
        {"=A,=A,=A": 10, "=A,+I,=A": 20},
        {"=C,=C,=C": 15, "=C,+I,=C": 25},
        {"=A,=A,=A": 10, "=A,+I,=A": 10},
        {"=A,=A,=A": 10, "=A,+I,=A": 40},
    ]
    mutation_loci = [set("+"), set(), set("+"), set("+")]
    result = score_handler.update_insertion_score(percent_discarded, mutation_loci)
    expected_result = [
        {"=A,=A,=A": 10, "=A,+I,=A": 20},
        {"=C,=C,=C": 15, "=C,+I,=C": 25},
        {"=A,=A,=A": 10, "=A,+I,=A": 40},
        {"=A,=A,=A": 10, "=A,+I,=A": 40},
    ]
    assert result == expected_result


###############################################################################
# annotate_score
###############################################################################
