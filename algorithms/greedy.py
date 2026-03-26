"""Weighted greedy spanner: sort by Euclidean distance, add edge iff w(p,q)*delta_E(p,q) > t*|pq|."""

from __future__ import annotations

from typing import Any, Optional

import numpy as np

from algorithms.base import (
    APSP_MODE_FULL,
    progress_iter,
    resolve_apsp_mode,
    resolve_candidates,
    resolve_weight,
)
from algorithms.registry import register
from core.metrics import (
    _build_adjacency,
    _euclidean_distances,
    _floyd_warshall,
    apsp_add_edge,
)


def _should_add_weighted_greedy(
    delta_E: float,
    w_pq: float,
    pq: float,
    t: float,
) -> bool:
    """Add edge iff w(p,q) * delta_E(p,q) > t * |pq|. When delta_E=inf (disconnected), add."""
    return np.isinf(delta_E) or (w_pq * delta_E > t * pq)


@register("greedy")
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
    Weighted greedy spanner algorithm.

    - Process candidate pairs sorted by Euclidean distance (ascending).
    - Add edge (p,q) iff w(p,q) * delta_E(p,q) > t * |pq|,
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

    # Sort by Euclidean distance (ascending)
    lengths = dist[candidates[:, 0], candidates[:, 1]].astype(np.float64)
    order = np.argsort(lengths)
    candidates = candidates[order]

    apsp_mode = resolve_apsp_mode(config)
    use_incremental = apsp_mode != APSP_MODE_FULL

    selected: list[tuple[int, int]] = []
    sp_dist: Optional[np.ndarray] = None  # APSP matrix for incremental mode

    for idx in progress_iter(
        range(candidates.shape[0]),
        total=candidates.shape[0],
        desc="greedy",
        config=config,
    ):
        i, j = int(candidates[idx, 0]), int(candidates[idx, 1])
        pq = dist[i, j]
        w_pq = w[i, j]

        if not selected:
            delta_E = np.inf
        elif use_incremental:
            delta_E = sp_dist[i, j]
        else:
            E_curr = np.array(selected, dtype=np.int64)
            adj = _build_adjacency(E_curr, n, dist)
            sp_dist = _floyd_warshall(adj)
            delta_E = sp_dist[i, j]

        if _should_add_weighted_greedy(delta_E, w_pq, pq, t):
            selected.append((i, j))
            if use_incremental:
                if sp_dist is None:
                    sp_dist = np.full((n, n), np.inf, dtype=np.float64)
                    np.fill_diagonal(sp_dist, 0.0)
                apsp_add_edge(sp_dist, i, j, pq)

    return np.array(selected, dtype=np.int64) if selected else np.empty((0, 2), dtype=np.int64)
