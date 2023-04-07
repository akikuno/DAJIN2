from src.DAJIN2.core.clustering.clustering import reorder_labels


def test_reorder_labels():
    labels = [1, 4, 4, 4, 2, 2, 2]
    test = reorder_labels(labels)
    answer = [1, 2, 2, 2, 3, 3, 3]
    assert test == answer


def test_reorder_labels_start_1():
    labels = [1, 4, 4, 4, 2, 2, 2]
    test = reorder_labels(labels, start=1)
    answer = [2, 3, 3, 3, 4, 4, 4]
    assert test == answer
