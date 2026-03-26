"""MST baseline weight for point sets."""

from __future__ import annotations

from typing import Optional

import numpy as np


def _euclidean_distances(points: np.ndarray) -> np.ndarray:
    """Compute pairwise Euclidean distances for points (n, 2). Returns n x n matrix."""
    n = points.shape[0]
    if n == 0:
        return np.empty((0, 0), dtype=np.float64)
    diff = points[:, np.newaxis, :] - points[np.newaxis, :, :]
    return np.sqrt(np.sum(diff**2, axis=2))


def _prim_mst(dist_matrix: np.ndarray) -> float:
    """
    Compute MST total weight via Prim's algorithm.
    dist_matrix is n x n symmetric; diagonal ignored.
    Returns 0.0 for n <= 1.
    """
    n = dist_matrix.shape[0]
    if n <= 1:
        return 0.0

    # Copy and ensure no self-edges
    D = np.array(dist_matrix, dtype=np.float64)
    np.fill_diagonal(D, np.inf)

    in_mst = np.zeros(n, dtype=bool)
    key = np.full(n, np.inf, dtype=np.float64)
    key[0] = 0.0
    total = 0.0

    for _ in range(n):
        u = np.argmin(key)
        if np.isinf(key[u]):
            break
        in_mst[u] = True
        total += key[u]
        key[u] = np.inf

        for v in range(n):
            if not in_mst[v] and D[u, v] < key[v]:
                key[v] = D[u, v]

    return total


def mst_weight(
    points: np.ndarray,
    weight: Optional[np.ndarray] = None,
    problem_type: str = "t",
) -> float:
    """
    Compute MST weight for point set (full graph on nodes).

    Args:
        points: Point coordinates, shape (n, 2).
        weight: n x n weight matrix for tw mode; None for t mode.
        problem_type: "t" (Euclidean distances) or "tw" (weight matrix).

    Returns:
        Total weight of minimum spanning tree. 0.0 for n <= 1.
    """
    n = points.shape[0]
    if n <= 1:
        return 0.0

    if problem_type == "t":
        dist_matrix = _euclidean_distances(points)
    else:
        if weight is None:
            weight = np.ones((n, n), dtype=np.float64)
        dist_matrix = np.asarray(weight, dtype=np.float64)

    return _prim_mst(dist_matrix)
