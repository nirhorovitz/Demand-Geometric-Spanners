"""
Weighted assignment greedy: propagation + refinement.

Intent: Propagation — spread weights along skeleton edges (symmetric, bounded).
       Refinement — remaining candidates ordered by m(p,q)=|pq|/w_prop(p,q),
       add iff delta_E > t*|pq|.

Stage 1: Greedy skeleton (Euclidean order).
Stage 2: Weight propagation along skeleton (config["propagation_delta"], default 0.1).
Stage 3: Weighted refinement over remaining candidates.
"""

from __future__ import annotations

from typing import Any, Optional

import numpy as np

from algorithms.base import progress_iter, resolve_candidates, resolve_weight
from algorithms.registry import register
from core.metrics import _build_adjacency, _euclidean_distances, _floyd_warshall

DEFAULT_PROPAGATION_DELTA = 0.1


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


def _propagate_weights(
    w: np.ndarray,
    skeleton: list[tuple[int, int]],
    n: int,
    delta: float,
    config: Optional[dict[str, Any]] = None,
) -> np.ndarray:
    """
    Propagate weights along skeleton edges. Symmetric and bounded.

    For each skeleton edge (a,b), add delta * w[a,b] to incident edges (a,k) and (b,k).
    Updates are symmetric (w[i,j] = w[j,i]) and bounded to [0, 1].

    Args:
        w: n x n weight matrix (input).
        skeleton: List of (i, j) edges in skeleton.
        n: Number of vertices.
        delta: Propagation strength in [0, 1].

    Returns:
        n x n propagated weight matrix.
    """
    w_prop = np.array(w, dtype=np.float64, copy=True)
    delta = max(0.0, min(1.0, float(delta)))

    for a, b in progress_iter(
        skeleton,
        total=len(skeleton),
        desc="assignment_greedy:propagate",
        config=config,
    ):
        wab = w_prop[a, b]
        for k in range(n):
            if k != a:
                val = w_prop[a, k] + delta * wab
                val = max(0.0, min(1.0, val))
                w_prop[a, k] = w_prop[k, a] = val
            if k != b:
                val = w_prop[b, k] + delta * wab
                val = max(0.0, min(1.0, val))
                w_prop[b, k] = w_prop[k, b] = val

    return w_prop


def _greedy_skeleton(
    candidates: np.ndarray,
    dist: np.ndarray,
    n: int,
    t: float,
    config: Optional[dict[str, Any]] = None,
) -> list[tuple[int, int]]:
    """
    Build greedy t-spanner skeleton ordered by Euclidean distance.
    Add edge (p,q) iff delta_E(p,q) > t*|pq|.
    """
    lengths = np.array([dist[int(i), int(j)] for i, j in candidates], dtype=np.float64)
    order = np.argsort(lengths)
    ordered = candidates[order]

    skeleton: list[tuple[int, int]] = []
    for idx in progress_iter(
        range(ordered.shape[0]),
        total=ordered.shape[0],
        desc="assignment_greedy:skeleton",
        config=config,
    ):
        i, j = int(ordered[idx, 0]), int(ordered[idx, 1])
        pq = dist[i, j]
        delta_E = np.inf
        if skeleton:
            E_curr = np.array(skeleton, dtype=np.int64)
            sp = _all_pairs_shortest_paths(E_curr, n, dist)
            delta_E = sp[i, j]
        if np.isinf(delta_E) or delta_E > t * pq:
            skeleton.append((i, j))
    return skeleton


@register("assignment_greedy")
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
    Weighted assignment greedy: propagation + refinement.

    Stage 1 (skeleton): Greedy t-spanner ordered by Euclidean distance.

    Stage 2 (propagation): Propagate weights along skeleton edges.
        Symmetric and bounded to [0,1]. config["propagation_delta"] (default 0.1).

    Stage 3 (refinement): Remaining candidates ordered by m(p,q)=|pq|/w_prop(p,q).
        Add edge iff delta_E > t*|pq|.

    E_input=None => full graph. weight=None => ones. Protocol preserved.
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

    cfg = config or {}
    delta = cfg.get("propagation_delta", DEFAULT_PROPAGATION_DELTA)
    delta = max(0.0, min(1.0, float(delta)))

    # Stage 1: Greedy skeleton
    skeleton = _greedy_skeleton(candidates, dist, n, t, config=config)
    skeleton_set = {tuple(sorted(e)) for e in skeleton}

    # Stage 2: Weight propagation (symmetric, bounded)
    w_prop = _propagate_weights(w, skeleton, n, delta, config=config)

    # Stage 3: Weighted refinement — remaining candidates ordered by m(p,q)
    remaining = [
        (int(c[0]), int(c[1]))
        for c in candidates
        if tuple(sorted((int(c[0]), int(c[1])))) not in skeleton_set
    ]
    if not remaining:
        return np.array(skeleton, dtype=np.int64)

    m_values = np.array(
        [
            dist[i, j] / max(w_prop[i, j], 1e-12)
            for i, j in remaining
        ],
        dtype=np.float64,
    )
    order_s3 = np.argsort(m_values)
    remaining_arr = np.array(remaining, dtype=np.int64)
    remaining_ordered = remaining_arr[order_s3]

    selected = list(skeleton)
    for idx in progress_iter(
        range(remaining_ordered.shape[0]),
        total=remaining_ordered.shape[0],
        desc="assignment_greedy:refine",
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
