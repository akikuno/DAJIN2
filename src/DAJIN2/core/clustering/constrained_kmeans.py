"""Constrained k-means clustering
Bradley, Paul S., Kristin P. Bennett, and Ayhan Demiriz. "Constrained k-means clustering." Microsoft Research, Redmond 20.0 (2000): 0.

The implementation is modified from @kuga-qiita. (https://qiita.com/kuga-qiita/items/5588d5469f3268b7fd39#constrained-k-means-1)
"""

from __future__ import annotations

from functools import partial

import numpy as np
from scipy.spatial.distance import cdist


def random(
    X: np.ndarray, n_clusters: int, random_state: np.random.RandomState, **kwargs
) -> tuple[np.ndarray, np.ndarray]:
    n_samples, _ = X.shape
    indices = random_state.choice(n_samples, size=n_clusters, replace=True)
    centers = X[indices]
    return centers, indices


def kmeans_plusplus(X, n_clusters, random_state):
    n_samples, _ = X.shape
    centers = np.empty((n_clusters, X.shape[1]), dtype=X.dtype)
    indices = np.empty(n_clusters, dtype=int)

    center_id = random_state.choice(n_samples)
    centers[0] = X[center_id]
    indices[0] = center_id

    closest_dist_sq = np.full(n_samples, np.inf)
    for i in range(1, n_clusters):
        dist_sq = np.sum((X - centers[i - 1]) ** 2, axis=1)
        closest_dist_sq = np.minimum(closest_dist_sq, dist_sq)
        probs = closest_dist_sq / np.sum(closest_dist_sq)
        cumulative_probs = np.cumsum(probs)
        r = random_state.rand()
        center_id = np.searchsorted(cumulative_probs, r)
        centers[i] = X[center_id]
        indices[i] = center_id

    return centers, indices


class ConstrainedKMeans:
    def __init__(
        self,
        n_clusters: int,
        min_cluster_size: int,
        max_cluster_size: int | None = None,
        random_state: int = 0,
        init: str = "k-means++",
        n_init: int = 10,
        max_iter: int = 300,
        tol: float = 1e-4,
        **kwargs,
    ):
        self.valid(
            n_clusters,
            min_cluster_size,
            max_cluster_size,
            random_state,
            init,
            n_init,
            max_iter,
            tol,
        )
        self.n_clusters = n_clusters
        self.min_cluster_size = min_cluster_size
        self.max_cluster_size = max_cluster_size
        self.random_state = np.random.RandomState(random_state)

        if init == "random":
            self.init = partial(
                random,
                n_clusters=self.n_clusters,
                random_state=self.random_state,
                **kwargs,
            )
        elif init == "k-means++":
            self.init = partial(
                kmeans_plusplus,
                n_clusters=self.n_clusters,
                random_state=self.random_state,
                **kwargs,
            )

        self.n_init = n_init
        self.max_iter = max_iter
        self.tol = tol
        self.kwargs = kwargs
        self.centers = None
        self.labels = None

    def get_lp_params(self, n_samples: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[float]]:
        X_nodes = np.arange(n_samples)
        cluster_nodes = np.arange(n_samples, n_samples + self.n_clusters)
        artificial_demand_node = np.array([n_samples + self.n_clusters])
        start_nodes = np.concatenate([np.repeat(X_nodes, self.n_clusters), cluster_nodes])
        end_nodes = np.concatenate(
            [
                np.tile(cluster_nodes, n_samples),
                np.repeat(artificial_demand_node, self.n_clusters),
            ]
        )
        capacities = np.concatenate(
            [
                np.ones(self.n_clusters * n_samples),
                np.full(
                    self.n_clusters,
                    n_samples - self.n_clusters * self.min_cluster_size
                    if self.max_cluster_size is None
                    else self.max_cluster_size - self.min_cluster_size,
                ),
            ]
        )
        supplies = np.concatenate(
            [
                np.ones(n_samples),
                -1 * np.full(self.n_clusters, self.min_cluster_size),
                np.array([-n_samples + self.n_clusters * self.min_cluster_size]),
            ]
        ).tolist()
        return start_nodes, end_nodes, capacities, supplies

    def simulated_annealing(self, X: np.ndarray, centers: np.ndarray) -> np.ndarray:
        """
        Perform simulated annealing to optimize cluster assignments.

        This method uses simulated annealing to assign data points to clusters,
        allowing for probabilistic acceptance of worse assignments to escape local optima.

        Parameters:
        X (np.ndarray): The input data array of shape (n_samples, n_features).
        centers (np.ndarray): The initial cluster centers of shape (n_clusters, n_features).

        Returns:
        np.ndarray: A mask array of shape (n_samples, n_clusters) indicating cluster assignments.
        """
        n_samples, n_clusters = X.shape[0], self.n_clusters
        dist_sq = cdist(X, centers, "sqeuclidean")
        labels = np.argmin(dist_sq, axis=1)

        T = 1.0
        T_min = 0.0001
        alpha = 0.9

        while T > T_min:
            # Vectorize sample processing and compute in batches
            idx = np.random.randint(0, n_samples, size=1000)
            current_labels = labels[idx]
            new_labels = np.random.randint(0, n_clusters, size=1000)

            # Select only the samples where the new label differs from the current label
            mask = current_labels != new_labels
            idx = idx[mask]
            current_labels = current_labels[mask]
            new_labels = new_labels[mask]

            if len(idx) == 0:
                T *= alpha
                continue

            current_dists = dist_sq[idx, current_labels]
            new_dists = dist_sq[idx, new_labels]

            delta_E = new_dists - current_dists
            prob = np.exp(-delta_E / T)
            accept = (delta_E < 0) | (prob > np.random.rand(len(prob)))

            labels[idx[accept]] = new_labels[accept]

            T *= alpha

        mask = np.zeros((n_samples, n_clusters))
        mask[np.arange(n_samples), labels] = 1
        return mask

    def fit_once(self, X: np.ndarray, **kwargs) -> tuple[np.ndarray, float]:
        centers_, _ = self.init(X, **kwargs)
        dist_sq = cdist(X, centers_, "sqeuclidean")
        mask = self.simulated_annealing(X, centers_)
        inertia = np.sum(mask * dist_sq)
        return centers_, inertia

    def fit(self, X: np.ndarray) -> ConstrainedKMeans:
        self.valid_tr(X)
        if X.ndim == 1:
            X = X[:, np.newaxis]

        best_inertia, best_centers = None, None
        for _ in range(self.n_init):
            centers_, inertia = self.fit_once(X, **self.kwargs)
            if best_inertia is None or inertia < best_inertia:
                best_inertia = inertia
                best_centers = centers_
        self.centers = best_centers

        dist_sq = cdist(X, self.centers, "sqeuclidean")

        for i in range(self.max_iter):
            mask = self.simulated_annealing(X, self.centers)
            centers_ = np.dot(mask.T, X) / np.sum(mask, axis=0)[:, np.newaxis]
            dist_sq = cdist(X, centers_, "sqeuclidean")
            centers_squared_diff = np.sum((centers_ - self.centers) ** 2)
            self.centers = centers_
            iter_ = i
            if centers_squared_diff <= self.tol:
                break

        self.inertia = np.sum(mask * dist_sq)
        self.labels = np.argmax(mask, axis=-1)
        self.iter_ = iter_
        return self

    def fit_predict(self, X: np.ndarray) -> np.ndarray:
        return self.fit(X).labels

    def predict(self, X: np.ndarray) -> np.ndarray:
        self.valid_pr(X)
        if X.ndim == 1:
            X = X[:, np.newaxis]
        mask = self.simulated_annealing(X, self.centers)
        labels = np.argmax(mask, axis=1)
        return labels

    def valid(
        self,
        n_clusters: int,
        min_cluster_size: int,
        max_cluster_size: int | None,
        random_state: int,
        init: str,
        n_init: int,
        max_iter: int,
        tol: float,
    ) -> None:
        if not isinstance(n_clusters, int):
            raise TypeError("n_clusters must be an integer")
        if not isinstance(min_cluster_size, int):
            raise TypeError("min_cluster_size must be an integer")
        if not isinstance(max_cluster_size, int) and max_cluster_size is not None:
            raise TypeError("max_cluster_size must be an integer or None")
        if not isinstance(random_state, int) and random_state is not None:
            raise TypeError("random_state must be an integer")
        if init not in ["random", "k-means++"]:
            raise TypeError("init must be 'random' or 'k-means++'")
        if not isinstance(n_init, int):
            raise TypeError("n_init must be an integer")
        if not isinstance(max_iter, int):
            raise TypeError("max_iter must be an integer")
        if not isinstance(tol, (int, float)):
            raise TypeError("tol must be an integer or float")

        if n_clusters < 2:
            raise ValueError("n_clusters must be 2 or more")
        if min_cluster_size < 0:
            raise ValueError("min_cluster_size must be 0 or more")
        if isinstance(max_cluster_size, int):
            if max_cluster_size < min_cluster_size or max_cluster_size < 0:
                raise ValueError("max_cluster_size must be at least min_cluster_size")
        if isinstance(random_state, int):
            if random_state < 0:
                raise ValueError("random_state must be 0 or more")
        if n_init <= 0:
            raise ValueError("n_init must be 1 or more")
        if max_iter < 1:
            raise ValueError("max_iter must be 1 or more")
        if tol < 0:
            raise ValueError("tol must be 0 or more")

    def valid_tr(self, X: np.ndarray) -> None:
        if not isinstance(X, np.ndarray):
            raise TypeError("X must be numpy.ndarray")
        if X.ndim < 1 or X.ndim > 2:
            raise ValueError("X must be 1-D or 2-D")
        if X.shape[0] < self.n_clusters or X.shape[0] < self.n_clusters * self.min_cluster_size:
            raise ValueError("X must have at least n_clusters * min_cluster_size rows")
        if self.max_cluster_size is not None:
            if X.shape[0] > self.n_clusters * self.max_cluster_size:
                raise ValueError("X must have at most n_clusters * max_cluster_size rows")

    def valid_pr(self, X: np.ndarray) -> None:
        if self.centers is None or self.labels is None:
            raise RuntimeError("fit method hasn't been called")
        if not isinstance(X, np.ndarray):
            raise TypeError("X must be numpy.ndarray")
        if X.ndim < 1 or X.ndim > 2:
            raise ValueError("X must be 1-D or 2-D")
        if X.shape[0] < self.n_clusters or X.shape[0] < self.n_clusters * self.min_cluster_size:
            raise ValueError("X must have at least n_clusters * min_cluster_size rows")
        if self.max_cluster_size is not None:
            if X.shape[0] > self.n_clusters * self.max_cluster_size:
                raise ValueError("X must have at most n_clusters * max_cluster_size rows")
        if X.ndim == 2 and X.shape[-1] != self.centers.shape[-1]:
            raise ValueError("input shape does not match the centers shape")
