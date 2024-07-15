"""Constrained k-means clustering
Bradley, Paul S., Kristin P. Bennett, and Ayhan Demiriz. "Constrained k-means clustering." Microsoft Research, Redmond 20.0 (2000): 0.

The implementation is referenced from @kuga-qiita. (https://qiita.com/kuga-qiita/items/5588d5469f3268b7fd39#constrained-k-means-1)
"""

from __future__ import annotations

from functools import partial

import numpy as np
from ortools.graph.python import min_cost_flow
from scipy.spatial.distance import cdist


def random(
    X: np.ndarray, n_clusters: int, random_state: np.random.RandomState, **kwargs
) -> tuple[np.ndarray, np.ndarray]:
    n_samples, _ = X.shape
    indices = random_state.choice(n_samples, size=n_clusters, replace=True)
    centers = X[indices]
    return centers, indices


class ConstrainedKMeans:
    def __init__(
        self,
        n_clusters: int,
        min_cluster_size: int,
        max_cluster_size: int | None = None,
        random_state: int = 0,
        init: str = "random",
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

        self.n_init = n_init
        self.max_iter = max_iter
        self.tol = tol
        self.kwargs = kwargs
        self.centers = None
        self.labels = None

    def get_smcf_params(self, n_samples: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[int]]:
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

    def calc_unit_costs(self, X: np.ndarray, centers: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        dist_sq = cdist(X, centers, "sqeuclidean")
        unit_costs = np.concatenate([dist_sq.flatten(), np.zeros(self.n_clusters)])
        return unit_costs, dist_sq

    def clustering(
        self,
        start_nodes: np.ndarray,
        end_nodes: np.ndarray,
        capacities: np.ndarray,
        supplies: list[int],
        unit_costs: np.ndarray,
    ) -> np.ndarray:
        smcf = min_cost_flow.SimpleMinCostFlow()
        all_arcs = smcf.add_arcs_with_capacity_and_unit_cost(start_nodes, end_nodes, capacities, unit_costs)
        smcf.set_nodes_supplies(np.arange(len(supplies)), supplies)
        smcf.solve()
        solution_flows = smcf.flows(all_arcs)
        mask = solution_flows[: -self.n_clusters].reshape((-1, self.n_clusters))
        return mask

    def fit(self, X: np.ndarray) -> ConstrainedKMeans:
        self.valid_tr(X)
        if X.ndim == 1:
            X = X[:, np.newaxis]
        n_samples, _ = X.shape
        params = self.get_smcf_params(n_samples)

        best_inertia, centers, self.init_indices = None, None, None
        for _ in range(self.n_init):
            centers_, init_indices_ = self.init(X, **self.kwargs)
            unit_costs, dist_sq = self.calc_unit_costs(X, centers_)
            mask = self.clustering(*params, unit_costs)
            inertia = np.sum(mask * dist_sq)
            if best_inertia is None or best_inertia > inertia:
                best_inertia = inertia
                centers = centers_
                self.init_indices = init_indices_

        unit_costs, dist_sq = self.calc_unit_costs(X, centers)

        for i in range(self.max_iter):
            mask = self.clustering(*params, unit_costs)
            centers_ = np.dot(mask.T, X) / np.sum(mask, axis=0)[:, np.newaxis]
            unit_costs, dist_sq = self.calc_unit_costs(X, centers_)
            centers_squared_diff = np.sum((centers_ - centers) ** 2)
            centers = centers_
            iter_ = i
            if centers_squared_diff <= self.tol:
                break

        self.inertia = np.sum(mask * dist_sq)
        self.centers = centers
        self.labels = np.argmax(mask, axis=-1)
        self.iter_ = iter_
        return self

    def fit_predict(self, X: np.ndarray) -> np.ndarray:
        return self.fit(X).labels

    def predict(self, X: np.ndarray) -> np.ndarray:
        self.valid_pr(X)
        if X.ndim == 1:
            X = X[:, np.newaxis]
        n_samples, _ = X.shape
        params = self.get_smcf_params(n_samples)
        unit_costs = self.calc_unit_costs(X, self.centers)
        mask = self.clustering(*params, unit_costs)
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
        if init not in ["random"]:
            raise TypeError("init must be 'random'")
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
