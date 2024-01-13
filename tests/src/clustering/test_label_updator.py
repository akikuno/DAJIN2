from src.DAJIN2.core.clustering import label_updator


def test_relabel_with_consective_order():
    labels = [1, 4, 4, 4, 2, 2, 2]
    test = label_updator.relabel_with_consective_order(labels)
    answer = [1, 2, 2, 2, 3, 3, 3]
    assert test == answer


def test_relabel_with_consective_order_start_2():
    labels = [1, 4, 4, 4, 2, 2, 2]
    test = label_updator.relabel_with_consective_order(labels, start=2)
    answer = [2, 3, 3, 3, 4, 4, 4]
    assert test == answer
