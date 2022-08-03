import scipy.stats as st
import numpy as np


def test_fisher_exact():
    data = [[4, 14], [3, 1]]
    test = st.fisher_exact(data)[1]
    answer = 0.07655502392344499
    assert np.isclose(test, answer)


def test_chisquare():
    data = [[4, 14], [3, 1]]
    test = st.chi2_contingency(data)[1]
    answer = 0.1452510053118262
    assert np.isclose(test, answer)
