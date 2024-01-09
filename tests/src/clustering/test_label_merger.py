from DAJIN2.core.clustering.label_merger import merge_mixed_cluster
from DAJIN2.core.clustering.label_merger import map_clusters_to_previous
from DAJIN2.core.clustering.label_merger import merge_minor_cluster


def test_merge_mixed_cluster():
    labels_control = [1, 1, 2, 3, 3, 3, 4, 4, 4]
    labels_sample = [1, 2, 2, 2, 3, 4, 5, 5]
    # label 1, 3, 4 are common in both control and sample so they should be merged
    expected_output = [6, 2, 2, 2, 6, 6, 5, 5]  # 6 is the new label for merged clusters
    assert merge_mixed_cluster(labels_control, labels_sample, threshold=20) == expected_output


def test_map_clusters_to_previous():
    labels_sample = [1, 1, 1, 0, 2, 0]
    labels_previous = [0, 0, 0, 1, 1, 1]
    expected_correspondence = {1: 0, 0: 1, 2: 1}

    assert map_clusters_to_previous(labels_sample, labels_previous) == expected_correspondence


def test_merge_minor_cluster():
    """
    In labels_sample, label 1 occurs 37.5% of the time, label 3 occurs 25%, and labels 2, 4, 5 occur 12.5% each. Therefore, labels 2, 4, 5 should be merged into previous lables
    """
    labels_sample = [1, 1, 1, 2, 3, 3, 4, 5]
    labels_previous = [0, 0, 0, 0, 1, 1, 2, 2]
    expected_output = [0, 0, 0, 0, 1, 1, 2, 2]
    assert (
        merge_minor_cluster(labels_sample, labels_previous, threshold_percentage=20, threshold_readnumber=100)
        == expected_output
    )
