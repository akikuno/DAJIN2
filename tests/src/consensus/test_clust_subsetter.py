from DAJIN2.core.consensus.clust_formatter import subset_clust


def test_basic_subset_clust():
    sample = [{"LABEL": 1, "DATA": i} for i in range(1500)]
    subset = subset_clust(sample)
    assert len(subset) == 1000  # Default num value is 1000
    assert all(item["LABEL"] == 1 for item in subset)


def test_subset_clust_with_custom_num():
    sample = [{"LABEL": 1, "DATA": i} for i in range(1500)]
    subset = subset_clust(sample, num=500)
    assert len(subset) == 500


def test_subset_clust_with_multiple_labels():
    sample = [{"LABEL": 1, "DATA": i} for i in range(500)] + [{"LABEL": 2, "DATA": i} for i in range(500)]
    subset = subset_clust(sample)
    assert len(subset) == 1000
    assert len([item for item in subset if item["LABEL"] == 1]) == 500
    assert len([item for item in subset if item["LABEL"] == 2]) == 500


def test_subset_clust_with_less_than_num():
    sample = [{"LABEL": 1, "DATA": i} for i in range(400)]
    subset = subset_clust(sample, num=500)
    assert len(subset) == 400


def test_subset_clust_with_multiple_labels_with_less_than_nu():
    sample = [{"LABEL": 1, "DATA": i} for i in range(500)] + [{"LABEL": 2, "DATA": i} for i in range(500)]
    subset = subset_clust(sample, num=200)
    assert len(subset) == 400
    assert len([item for item in subset if item["LABEL"] == 1]) == 200
    assert len([item for item in subset if item["LABEL"] == 2]) == 200


def test_subset_clust_with_empty_sample():
    sample = []
    subset = subset_clust(sample)
    assert len(subset) == 0
