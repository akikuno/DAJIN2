import numpy as np
from sklearn.neighbors import LocalOutlierFactor as lof

X = [[-1.1], [0.2], [101.1], [0.3]]
clf = lof(n_neighbors=2)
clf.fit_predict(X)

clf.negative_outlier_factor_

X = [[-1.1], [0.2], [0.5], [0.3]]
clf = lof(n_neighbors=2, novelty=True)
clf.fit(X)
clf.predict([[0], [101.1]])
