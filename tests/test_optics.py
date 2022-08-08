from sklearn.cluster import Birch
from sklearn.cluster import OPTICS
import numpy as np

X = np.array([[1, 2], [2, 5], [3, 6], [8, 7], [8, 8], [7, 3]])
clustering = OPTICS(min_samples=2).fit(X)
clustering.labels_


X = [[100, 100], [0, 1], [0.3, 1], [-0.3, 1], [0, -1], [0.3, -1], [-0.3, -1]]
brc = Birch(n_clusters=None)
brc.fit(X)

brc.predict(X)
