from DAJIN2.core.consensus.clust_formatter import downsample_by_label


def test_basic_downsample_by_label():
    sample = [{"LABEL": 1, "DATA": i} for i in range(1500)]
    subset = downsample_by_label(sample)
    assert len(subset) == 1000  # Default num value is 1000
    assert all(item["LABEL"] == 1 for item in subset)


def test_downsample_by_label_with_custom_num():
    sample = [{"LABEL": 1, "DATA": i} for i in range(1500)]
    subset = downsample_by_label(sample, num=500)
    assert len(subset) == 500


def test_downsample_by_label_with_multiple_labels():
    sample = [{"LABEL": 1, "DATA": i} for i in range(500)] + [{"LABEL": 2, "DATA": i} for i in range(500)]
    subset = downsample_by_label(sample)
    assert len(subset) == 1000
    assert len([item for item in subset if item["LABEL"] == 1]) == 500
    assert len([item for item in subset if item["LABEL"] == 2]) == 500


def test_downsample_by_label_with_less_than_num():
    sample = [{"LABEL": 1, "DATA": i} for i in range(400)]
    subset = downsample_by_label(sample, num=500)
    assert len(subset) == 400


def test_downsample_by_label_with_multiple_labels_with_less_than_nu():
    sample = [{"LABEL": 1, "DATA": i} for i in range(500)] + [{"LABEL": 2, "DATA": i} for i in range(500)]
    subset = downsample_by_label(sample, num=200)
    assert len(subset) == 400
    assert len([item for item in subset if item["LABEL"] == 1]) == 200
    assert len([item for item in subset if item["LABEL"] == 2]) == 200


def test_downsample_by_label_with_empty_sample():
    sample = []
    subset = downsample_by_label(sample)
    assert len(subset) == 0
