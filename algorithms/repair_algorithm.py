"""
Repair-style iterative augmentation from greedy seed.

Intent: Iterative repair — greedy skeleton + repeated repair passes until
convergence or hard cap. Each pass considers remaining candidates ordered by
m(p,q)=|pq|/w(p,q), adding edges that improve stretch (delta_E > t*|pq|).

Runtime-safe loop controls (mandatory):
- max_iterations: config["max_iterations"] or min(100, 2*n)
- Stop immediately when a pass adds zero edges (convergence)
"""

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

# Termination safeguards: prevent non-terminating loops
DEFAULT_MAX_ITERATIONS = 100
DEFAULT_MAX_ITERATIONS_PER_N = 2  # max_iter = min(default, n * this)


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
    apsp_mode = resolve_apsp_mode(config)
    use_incremental = apsp_mode != APSP_MODE_FULL

    lengths = dist[candidates[:, 0], candidates[:, 1]].astype(np.float64)
    order = np.argsort(lengths)
    ordered = candidates[order]

    skeleton: list[tuple[int, int]] = []
    sp_dist: Optional[np.ndarray] = None

    for idx in progress_iter(
        range(ordered.shape[0]),
        total=ordered.shape[0],
        desc="repair:skeleton",
        config=config,
    ):
        i, j = int(ordered[idx, 0]), int(ordered[idx, 1])
        pq = dist[i, j]
        if not skeleton:
            delta_E = np.inf
        elif use_incremental:
            delta_E = sp_dist[i, j]
        else:
            E_curr = np.array(skeleton, dtype=np.int64)
            adj = _build_adjacency(E_curr, n, dist)
            sp_dist = _floyd_warshall(adj)
            delta_E = sp_dist[i, j]
        if np.isinf(delta_E) or delta_E > t * pq:
            skeleton.append((i, j))
            if use_incremental:
                if sp_dist is None:
                    sp_dist = np.full((n, n), np.inf, dtype=np.float64)
                    np.fill_diagonal(sp_dist, 0.0)
                apsp_add_edge(sp_dist, i, j, pq)
    return skeleton


def _repair_pass(
    selected: list[tuple[int, int]],
    remaining: list[tuple[int, int]],
    dist: np.ndarray,
    w: np.ndarray,
    n: int,
    t: float,
    config: Optional[dict[str, Any]] = None,
) -> tuple[list[tuple[int, int]], int]:
    """
    One repair pass: consider remaining candidates ordered by m(p,q)=|pq|/w(p,q),
    add edge iff delta_E > t*|pq|. Returns (updated selected, count of edges added).
    """
    if not remaining:
        return selected, 0

    apsp_mode = resolve_apsp_mode(config)
    use_incremental = apsp_mode != APSP_MODE_FULL

    rem_arr = np.array(remaining, dtype=np.int64)
    w_vals = np.maximum(w[rem_arr[:, 0], rem_arr[:, 1]], 1e-12)
    m_values = (dist[rem_arr[:, 0], rem_arr[:, 1]] / w_vals).astype(np.float64)
    order = np.argsort(m_values)
    ordered = rem_arr[order]

    added = 0
    sp_dist: Optional[np.ndarray] = None
    if use_incremental and selected:
        sp_dist = _build_adjacency(np.array(selected, dtype=np.int64), n, dist)
        sp_dist = _floyd_warshall(sp_dist)

    for idx in progress_iter(
        range(ordered.shape[0]),
        total=ordered.shape[0],
        desc="repair:pass",
        config=config,
    ):
        i, j = int(ordered[idx, 0]), int(ordered[idx, 1])
        pq = dist[i, j]
        if not selected:
            delta_E = np.inf
        elif use_incremental and sp_dist is not None:
            delta_E = sp_dist[i, j]
        else:
            E_curr = np.array(selected, dtype=np.int64)
            adj = _build_adjacency(E_curr, n, dist)
            sp_dist = _floyd_warshall(adj)
            delta_E = sp_dist[i, j]
        if np.isinf(delta_E) or delta_E > t * pq:
            selected.append((i, j))
            added += 1
            if use_incremental:
                if sp_dist is None:
                    sp_dist = np.full((n, n), np.inf, dtype=np.float64)
                    np.fill_diagonal(sp_dist, 0.0)
                apsp_add_edge(sp_dist, i, j, pq)
    return selected, added


@register("repair")
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
    Repair-style iterative augmentation from greedy seed.

    Stage 1 (seed): Greedy t-spanner ordered by Euclidean distance.
        Insertion: add edge (p,q) iff delta_E(p,q) > t*|pq|.

    Stage 2 (iterative repair): Consider remaining candidates ordered by
        m(p,q)=|pq|/w(p,q). Add edge iff delta_E > t*|pq|. Repeat until
        no new edges added (convergence) or max_iterations reached.

    Runtime-safe loop controls:
    - max_iterations: config["max_iterations"] or min(100, 2*n)
    - Stop when pass adds zero edges (convergence)

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

    # Resolve max_iterations from config (termination safeguard)
    cfg = config or {}
    max_iter = cfg.get("max_iterations")
    if max_iter is None:
        max_iter = min(DEFAULT_MAX_ITERATIONS, max(1, n * DEFAULT_MAX_ITERATIONS_PER_N))
    max_iter = max(1, int(max_iter))
    print(f"max_iter: {max_iter}")

    # Stage 1: Greedy skeleton
    skeleton = _greedy_skeleton(candidates, dist, n, 2, config=config)
    skeleton_set = {tuple(sorted(e)) for e in skeleton}

    # Remaining candidates (not in skeleton)
    remaining = [
        (int(c[0]), int(c[1]))
        for c in candidates
        if tuple(sorted((int(c[0]), int(c[1])))) not in skeleton_set
    ]
    if not remaining:
        return np.array(skeleton, dtype=np.int64)

    # Stage 2: Iterative repair loop — runtime-safe: convergence or max_iter
    selected = list(skeleton)
    for _ in range(max_iter):
        _, added = _repair_pass(selected, remaining, dist, w, n, t, config=config)
        if added == 0:
            break  # Convergence: no improvement this pass
        selected_set = {tuple(sorted(e)) for e in selected}
        remaining = [
            (i, j) for i, j in remaining
            if tuple(sorted((i, j))) not in selected_set
        ]
        if not remaining:
            break

    return np.array(selected, dtype=np.int64) if selected else np.empty((0, 2), dtype=np.int64)
