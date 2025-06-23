"""
Anomaly detection utilities for mutation extraction.
"""

from __future__ import annotations

import numpy as np
from sklearn.neural_network import MLPClassifier


def cosine_distance(x: list[float], y: list[float]) -> float:
    """Calculate cosine distance between two vectors."""
    # Add 1e-6 to avoid division by zero when calculating cosine similarity
    x = np.array(x) + 1e-6
    y = np.array(y) + 1e-6
    return 1 - float(np.dot(x, y) / (np.linalg.norm(x) * np.linalg.norm(y)))


def is_dissimilar_loci(values_sample, values_control, index: int, is_consensus: bool = False) -> bool:
    """Check if a locus shows dissimilar mutation patterns between sample and control."""
    # If 'sample' has more than 20% variation compared to 'control' in consensus mode,
    # unconditionally set it to 'dissimilar loci'. This is set to counteract cases where,
    # when evaluating cosine similarity during significant deletions, values exceedingly
    # close to 1 can occur even if not observed in the control
    # (e.g., control = [1,1,1,1,1], sample = [100,100,100,100,100] -> cosine similarity = 1).
    if values_sample[index] - values_control[index] > 20:
        if not is_consensus or values_sample[index] > 50:
            return True
        else:
            return False

    # Subset 10 bases around index.
    x = values_sample[index : index + 10]
    y = values_control[index : index + 10]

    x_slice = x[1:]
    y_slice = y[1:]

    distance = cosine_distance(x, y)
    distance_slice = cosine_distance(x_slice, y_slice)

    return distance > 0.05 and distance / (distance + distance_slice) > 0.9


def detect_anomalies(values_sample, values_control, threshold: float, is_consensus: bool = False) -> set[int]:
    """
    Detect anomalies and return indices of outliers.
    """
    rng = np.random.default_rng(seed=1)

    random_size = 10_000
    control_size = len(values_control)
    total_size = random_size + control_size

    randoms = rng.uniform(0, 100, random_size)
    randoms_error = np.clip(randoms + rng.uniform(0, threshold, random_size), 0, 100)
    randoms_mutation = np.clip(randoms + rng.uniform(threshold, 100, random_size), 0, 100)

    values_error = np.clip(values_control + rng.uniform(0, threshold, control_size), 0, 100)
    values_mutation = np.clip(values_control + rng.uniform(threshold, 100, control_size), 0, 100)

    matrix_error_randoms = np.array([randoms, randoms_error]).T
    matrix_error_control = np.array([values_control, values_error]).T
    matrix_error = np.concatenate([matrix_error_randoms, matrix_error_control], axis=0)

    matrix_mutation_randoms = np.array([randoms, randoms_mutation]).T
    matrix_mutation_control = np.array([values_control, values_mutation]).T
    matrix_mutation = np.concatenate([matrix_mutation_randoms, matrix_mutation_control], axis=0)

    X = np.concatenate([matrix_error, matrix_mutation], axis=0)
    y = [0] * (total_size) + [1] * (total_size)

    clf = MLPClassifier(solver="lbfgs", alpha=1e-5, hidden_layer_sizes=(5, 2), random_state=1)
    clf.fit(X, y)

    results = clf.predict(np.array([values_control, values_sample]).T)

    candidate_loci = {i for i, v in enumerate(results) if v == 1}

    return {i for i in candidate_loci if is_dissimilar_loci(values_sample, values_control, i, is_consensus)}


def extract_anomal_loci(
    indels_normalized_sample,
    indels_normalized_control,
    thresholds: dict[str, float],
    is_consensus: bool = False,
) -> dict[str, set[int]]:
    """Extract outlier loci comparing indel counts between sample and control."""
    anomal_loci = {}
    for mut in {"+", "-", "*"}:
        values_sample = indels_normalized_sample[mut]
        values_control = indels_normalized_control[mut]
        idx_outliers = detect_anomalies(values_sample, values_control, thresholds[mut], is_consensus)
        anomal_loci[mut] = idx_outliers
    return anomal_loci
