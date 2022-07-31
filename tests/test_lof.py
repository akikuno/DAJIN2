import numpy as np
from sklearn.neighbors import LocalOutlierFactor as lof


def test_lof_example():
    X = [[-1.1], [0.2], [101.1], [0.3]]
    clf = lof(n_neighbors=2)
    clf.fit_predict(X)
    test = clf.negative_outlier_factor_
    answer = np.array([-0.98214286, -1.03703704, -73.36970899, -0.98214286])
    assert np.allclose(test, answer)


def test_lof_novelty_example():
    X = [[-1.1], [0.2], [0.5], [0.3]]
    clf = lof(n_neighbors=2, novelty=True)
    clf.fit(X)
    test = clf.predict([[0], [101.1]])
    answer = np.array([1, -1])
    assert np.allclose(test, answer)
