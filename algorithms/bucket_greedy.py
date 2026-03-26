"""Bucket greedy: partition candidates by distance buckets (alpha>1), sort within bucket by weight, apply insertion condition."""

from __future__ import annotations

import math
from typing import Any, Optional

import numpy as np

from algorithms.base import progress_iter, resolve_candidates, resolve_weight
from algorithms.registry import register
from core.metrics import _build_adjacency, _euclidean_distances, _floyd_warshall

DEFAULT_ALPHA = 2.0
MIN_ALPHA = 1.0 + 1e-9
_DIST_EPS = 1e-12
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


def _resolve_alpha(config: Optional[dict[str, Any]]) -> float:
    """Resolve alpha from config with safe default. alpha must be > 1."""
    if config is None:
        return DEFAULT_ALPHA
    alpha = config.get("alpha", DEFAULT_ALPHA)
    alpha = float(alpha)
    if alpha <= MIN_ALPHA:
        return DEFAULT_ALPHA
    return alpha


def _bucket_index(dist: float, alpha: float, eps: float = _DIST_EPS) -> int:
    """Compute bucket index for distance. Bucket k = [alpha^k, alpha^(k+1))."""
    d = max(dist, eps)
    return int(math.floor(math.log(d) / math.log(alpha)))


@register("bucket_greedy")
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
    Bucket greedy spanner algorithm.

    - Partition candidate pairs into distance buckets using alpha > 1.
      Bucket k contains edges with |pq| in [alpha^k, alpha^(k+1)).
    - Process buckets in ascending order (shortest distances first).
    - Within each bucket, sort by m(p,q) = |pq|/w(p,q) non-decreasing.
      Ties broken by (i,j) for determinism.
    - Insertion: add edge (p,q) iff delta_E(p,q) > t*|pq|.
    - config["alpha"] with safe fallback; alpha <= 1 => default 2.0.
    - E_input=None => full graph candidates. weight=None => ones matrix.
    """
    _ = rng_seed

    n = points.shape[0]
    if n <= 1:
        return np.empty((0, 2), dtype=np.int64)

    candidates = resolve_candidates(points, E_input)
    w = resolve_weight(weight, n)
    dist = _euclidean_distances(points)
    alpha = _resolve_alpha(config)

    if candidates.shape[0] == 0:
        return np.empty((0, 2), dtype=np.int64)

    lengths = np.array(
        [dist[int(i), int(j)] for i, j in candidates],
        dtype=np.float64,
    )
    m_values = np.array(
        [
            dist[int(i), int(j)] / max(w[int(i), int(j)], _WEIGHT_EPS)
            for i, j in candidates
        ],
        dtype=np.float64,
    )
    bucket_indices = np.array(
        [_bucket_index(d, alpha) for d in lengths],
        dtype=np.int64,
    )

    # Group by bucket
    bucket_to_indices: dict[int, list[int]] = {}
    for idx in progress_iter(
        range(len(candidates)),
        total=len(candidates),
        desc="bucket_greedy:bucketize",
        config=config,
    ):
        b = int(bucket_indices[idx])
        bucket_to_indices.setdefault(b, []).append(idx)

    # Process buckets in ascending order; within bucket sort by m, then index
    ordered_indices: list[int] = []
    for b in sorted(bucket_to_indices.keys()):
        indices = bucket_to_indices[b]
        indices_sorted = sorted(indices, key=lambda i: (m_values[i], i))
        ordered_indices.extend(indices_sorted)

    candidates_ordered = candidates[ordered_indices]

    selected: list[tuple[int, int]] = []
    for idx in progress_iter(
        range(candidates_ordered.shape[0]),
        total=candidates_ordered.shape[0],
        desc="bucket_greedy:select",
        config=config,
    ):
        i, j = int(candidates_ordered[idx, 0]), int(candidates_ordered[idx, 1])
        pq = dist[i, j]
        delta_E = np.inf
        if selected:
            E_curr = np.array(selected, dtype=np.int64)
            sp = _all_pairs_shortest_paths(E_curr, n, dist)
            delta_E = sp[i, j]

        if np.isinf(delta_E) or delta_E > t * pq:
            selected.append((i, j))

    return np.array(selected, dtype=np.int64) if selected else np.empty((0, 2), dtype=np.int64)
