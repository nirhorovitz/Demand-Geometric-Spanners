"""Double greedy: two-stage approach — greedy skeleton plus weighted refinement."""

from __future__ import annotations

from typing import Any, Optional

import numpy as np

from algorithms.base import progress_iter, resolve_candidates, resolve_weight
from algorithms.registry import register
from core.metrics import _build_adjacency, _euclidean_distances, _floyd_warshall

# Practical safeguard: avoid division by zero in m(p,q) = |pq|/w(p,q)
_WEIGHT_EPS = 1e-12


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


def _sort_by_distance_then_index(
    candidates: np.ndarray,
    lengths: np.ndarray,
) -> np.ndarray:
    """Sort by distance ascending; ties broken by (i, j) for determinism."""
    n_cand = candidates.shape[0]
    secondary = np.arange(n_cand, dtype=np.int64)
    order = np.lexsort((secondary, lengths))
    return candidates[order]


def _sort_by_m_then_index(
    remaining: np.ndarray,
    m_values: np.ndarray,
) -> np.ndarray:
    """Sort by m(p,q) ascending; ties broken by (i, j) for determinism."""
    n_rem = remaining.shape[0]
    secondary = np.arange(n_rem, dtype=np.int64)
    order = np.lexsort((secondary, m_values))
    return remaining[order]


@register("double_greedy")
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
    Two-stage double-greedy spanner algorithm.

    Stage 1 (skeleton): Build greedy t-spanner ordered by Euclidean distance.
        Insertion: add edge (p,q) iff delta_E(p,q) > t*|pq|.
        Ties broken by (i,j) for deterministic ordering.

    Stage 2 (refinement): Process remaining candidates ordered by
        m(p,q) = |pq|/w(p,q) non-decreasing. Add edge iff delta_E(p,q) > t*|pq|.
        Weights affect ordering only; same insertion condition.
        Ties broken by (i,j) for determinism.

    Practical safeguards:
    - weight < eps => use eps to avoid division by zero in m(p,q).
    - E_input=None => full graph candidates. weight=None => ones matrix.
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

    # --- Stage 1: Greedy skeleton (ordered by Euclidean distance) ---
    lengths = np.array(
        [dist[int(i), int(j)] for i, j in candidates],
        dtype=np.float64,
    )
    candidates_s1 = _sort_by_distance_then_index(candidates, lengths)

    skeleton: list[tuple[int, int]] = []
    for idx in progress_iter(
        range(candidates_s1.shape[0]),
        total=candidates_s1.shape[0],
        desc="double_greedy:skeleton",
        config=config,
    ):
        i, j = int(candidates_s1[idx, 0]), int(candidates_s1[idx, 1])
        pq = dist[i, j]
        delta_E = np.inf
        if skeleton:
            E_curr = np.array(skeleton, dtype=np.int64)
            sp = _all_pairs_shortest_paths(E_curr, n, dist)
            delta_E = sp[i, j]
        if np.isinf(delta_E) or delta_E > t * pq:
            skeleton.append((i, j))

    skeleton_set = {tuple(sorted(e)) for e in skeleton}

    # --- Stage 2: Weighted refinement ---
    remaining = [
        (int(c[0]), int(c[1]))
        for c in candidates
        if tuple(sorted((int(c[0]), int(c[1])))) not in skeleton_set
    ]
    if not remaining:
        return np.array(skeleton, dtype=np.int64)

    m_values = np.array(
        [
            dist[i, j] / max(w[i, j], _WEIGHT_EPS)
            for i, j in remaining
        ],
        dtype=np.float64,
    )
    remaining_arr = np.array(remaining, dtype=np.int64)
    remaining_ordered = _sort_by_m_then_index(remaining_arr, m_values)

    selected = list(skeleton)
    for idx in progress_iter(
        range(remaining_ordered.shape[0]),
        total=remaining_ordered.shape[0],
        desc="double_greedy:refine",
        config=config,
    ):
        i, j = int(remaining_ordered[idx, 0]), int(remaining_ordered[idx, 1])
        pq = dist[i, j]
        E_curr = np.array(selected, dtype=np.int64)
        sp = _all_pairs_shortest_paths(E_curr, n, dist)
        delta_E = sp[i, j]
        if np.isinf(delta_E) or delta_E > t * pq:
            selected.append((i, j))

    return np.array(selected, dtype=np.int64) if selected else np.empty((0, 2), dtype=np.int64)
