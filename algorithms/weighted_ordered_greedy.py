"""Weighted ordered greedy: sort by m(p,q)=|pq|/w(p,q) non-decreasing, add edge iff delta_E(p,q) > t*|pq|."""

from __future__ import annotations

from typing import Any, Optional

import numpy as np

from algorithms.base import progress_iter, resolve_candidates, resolve_weight
from algorithms.registry import register
from core.metrics import _build_adjacency, _euclidean_distances, _floyd_warshall


def _all_pairs_shortest_paths(
    edges: np.ndarray,
    n: int,
    dist_matrix: np.ndarray,
) -> np.ndarray:
    """All-pairs shortest paths in graph with given edge lengths. Returns n x n matrix."""
    if n == 0:
        return np.empty((0, 0), dtype=np.float64)
    adj = _build_adjacency(edges, n, dist_matrix)
    return _floyd_warshall(adj)


def _should_add_weighted_ordered(delta_E: float, pq: float, t: float) -> bool:
    """Add edge iff delta_E(p,q) > t * |pq|. When delta_E=inf (disconnected), add."""
    return np.isinf(delta_E) or (delta_E > t * pq)


@register("weighted_ordered_greedy")
def run(
    points: np.ndarray,
    t: float,
    *,
    E_input: Optional[np.ndarray] = None,
    weight: Optional[np.ndarray] = None,
    config: Optional[dict[str, Any]] = None,
    rng_seed: Optional[int] = None,
) -> np.ndarray:
    """
    Weighted ordered greedy spanner algorithm.

    - Process pairs sorted by m(p,q) = |pq| / w(p,q) non-decreasing.
    - Add edge (p,q) iff delta_E(p,q) > t * |pq| (normal greedy insertion),
      where delta_E(p,q) is shortest path in current graph from p to q (Euclidean lengths).
    - weight=None => ones matrix (resolve_weight).
    - E_input=None => full graph candidates.
    """
    _ = rng_seed

    n = points.shape[0]
    if n <= 1:
        return np.empty((0, 2), dtype=np.int64)

    candidates = resolve_candidates(points, E_input)
    w = resolve_weight(weight, n)
    dist = _euclidean_distances(points)

    if candidates.shape[0] == 0:
        return np.empty((0, 2), dtype=np.int64)

    # m(p,q) = |pq| / w(p,q), sort non-decreasing by m
    m_values = np.array(
        [
            dist[int(i), int(j)] / max(w[int(i), int(j)], 1e-12)
            for i, j in candidates
        ],
        dtype=np.float64,
    )
    order = np.argsort(m_values)
    candidates = candidates[order]

    selected: list[tuple[int, int]] = []
    for idx in progress_iter(
        range(candidates.shape[0]),
        total=candidates.shape[0],
        desc="weighted_ordered_greedy",
        config=config,
    ):
        i, j = int(candidates[idx, 0]), int(candidates[idx, 1])
        pq = dist[i, j]

        E_curr = np.array(selected, dtype=np.int64) if selected else np.empty((0, 2), dtype=np.int64)
        sp = _all_pairs_shortest_paths(E_curr, n, dist)
        delta_E = sp[i, j]

        if _should_add_weighted_ordered(delta_E, pq, t):
            selected.append((i, j))

    return np.array(selected, dtype=np.int64) if selected else np.empty((0, 2), dtype=np.int64)
