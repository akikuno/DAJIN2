from src.DAJIN2.core.clustering.merge_clusters import merge_mixed_cluster
from src.DAJIN2.core.clustering.merge_clusters import merge_minor_cluster


def test_merge_mixed_cluster():
    """
    Label 2 and 3 are common in both control and sample.
    """
    labels_control = [1] * 10 + [2] * 5 + [3] * 5
    labels_sample = [2] * 5 + [3] * 5 + [4] * 10
    test = merge_mixed_cluster(labels_control, labels_sample)
    answer = [5] * 10 + [4] * 10
    assert test == answer


# def test_merge_minor_cluster():
#     labels = [1] * 100
#     labels += [2, 3, 4]
#     labels += [5] * 100
#     test = merge_minor_cluster(labels)
#     answer = [1] * 100
#     answer += [6] * 3
#     answer += [5] * 100
#     assert test == answer

