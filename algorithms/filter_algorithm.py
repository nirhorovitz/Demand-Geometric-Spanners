"""
Filter algorithm: greedy seed + remove/reinsert validation.

Intent: Remove/reinsert validation — for each edge, attempt removal; keep
removed iff w(p,q)*delta_E' <= t*|pq| for all p,q \in P (graph without edge still valid).
Edges processed in descending order by length.
"""

from __future__ import annotations

import os
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
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
    apsp_shortest_path_after_removal,
)


# filter exists-check config (SP-50)
DEFAULT_FILTER_EXISTS_THREADED = True
DEFAULT_FILTER_EXISTS_MAX_WORKERS = 4
DEFAULT_FILTER_EXISTS_CHUNK_SIZE = 64


def all_pairs_for_validation_check(n: int) -> list[tuple[int, int]]:
    """
    Return unordered full-point-set pairs for validation check.
    Domain: all (r,s) with r < s over points [0..n-1]. Independent from E_input.
    """
    if n < 2:
        return []
    return [(r, s) for r in range(n) for s in range(r + 1, n)]


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


def _exists_violation_sequential(
    pairs: list[tuple[int, int]],
    sp: np.ndarray,
    dist: np.ndarray,
    w: np.ndarray,
    n: int,
    t: float,
) -> bool:
    """
    Check pairs sequentially; return True if any (r,s) violates delta(rs)*w(rs) > t*|rs|.
    Break-both: stop on first violation.
    """
    for r, s in pairs:
        delta_rs = sp[r, s]
        wij = w[r, s]
        rs_len = dist[r, s]
        if np.isinf(delta_rs):
            return True  # Disconnected => violation
        if wij * delta_rs > t * rs_len:
            return True
    return False


def _exists_violation_chunk(
    pair_indices: list[int],
    all_pairs: list[tuple[int, int]],
    sp: np.ndarray,
    dist: np.ndarray,
    w: np.ndarray,
    n: int,
    t: float,
    stop_signal: Optional[threading.Event] = None,
) -> bool:
    """
    Check a chunk of pairs; return True if any violation found.
    When stop_signal is set, bail early (best-effort cancellation).
    """
    for idx in pair_indices:
        if stop_signal is not None and stop_signal.is_set():
            return False  # Early exit; no violation found in this chunk
        r, s = all_pairs[idx]
        delta_rs = sp[r, s]
        wij = w[r, s]
        rs_len = dist[r, s]
        if np.isinf(delta_rs):
            return True
        if wij * delta_rs > t * rs_len:
            return True
    return False


def _run_threaded_exists_check(
    pairs: list[tuple[int, int]],
    sp: np.ndarray,
    dist: np.ndarray,
    w: np.ndarray,
    n: int,
    t: float,
    max_workers: int,
    chunk_size: int,
) -> bool:
    """Run threaded exists-check; return True if any violation found."""
    chunks = []
    for start in range(0, len(pairs), chunk_size):
        end = min(start + chunk_size, len(pairs))
        chunks.append(list(range(start, end)))
    executor = ThreadPoolExecutor(max_workers=max_workers)
    try:
        futures = {
            executor.submit(
                _exists_violation_chunk,
                chunk,
                pairs,
                sp,
                dist,
                w,
                n,
                t,
            ): chunk
            for chunk in chunks
        }
        for future in as_completed(futures):
            if future.result():
                return True
        return False
    finally:
        executor.shutdown(wait=True)


def _run_threaded_exists_check_break_both(
    pairs: list[tuple[int, int]],
    sp: np.ndarray,
    dist: np.ndarray,
    w: np.ndarray,
    n: int,
    t: float,
    max_workers: int,
    chunk_size: int,
) -> bool:
    """
    Run threaded exists-check.
    Stop as soon as first violation is found; best-effort cancellation of remaining work.
    Deterministic pair ordering via chunk order; semantic equivalence with sequential path.
    """
    chunks = []
    for start in range(0, len(pairs), chunk_size):
        end = min(start + chunk_size, len(pairs))
        chunks.append(list(range(start, end)))
    stop_signal = threading.Event()
    executor = ThreadPoolExecutor(max_workers=max_workers)
    try:
        futures = [
            executor.submit(
                _exists_violation_chunk,
                chunk,
                pairs,
                sp,
                dist,
                w,
                n,
                t,
                stop_signal,
            )
            for chunk in chunks
        ]
        for future in futures:
            if stop_signal.is_set():
                future.cancel()
                continue
            try:
                if future.result():
                    stop_signal.set()
                    for f in futures:
                        if f != future:
                            f.cancel()
                    return True
            except Exception:
                raise  # Propagate to caller for sequential fallback
        return False
    finally:
        executor.shutdown(wait=True)


def _exists_violation_pairs(
    without: list[tuple[int, int]],
    dist: np.ndarray,
    w: np.ndarray,
    n: int,
    t: float,
    config: Optional[dict[str, Any]],
) -> bool:
    """
    Return True iff any unordered pair (r,s) with r<s violates delta(rs)*w(rs) > t*|rs|
    in graph `without`.
    Independent from E_input subset semantics.
    Break-both: stop scan on first violation.
    Uses ThreadPoolExecutor when filter_exists_threaded=True; sequential fallback on error.
    Config: filter_exists_threaded (default True), filter_exists_max_workers,
    filter_exists_chunk_size, filter_exists_break_both (default True, SP-58).
    """
    if n < 2:
        return False
    pairs = all_pairs_for_validation_check(n)
    if not pairs:
        return False

    E_without = np.array(without, dtype=np.int64)
    sp = _all_pairs_shortest_paths(E_without, n, dist)

    cfg = config or {}
    use_threaded = cfg.get("filter_exists_threaded", DEFAULT_FILTER_EXISTS_THREADED)
    max_workers = cfg.get("filter_exists_max_workers")
    if max_workers is None:
        try:
            max_workers = len(os.sched_getaffinity(0))
        except AttributeError:
            max_workers = os.cpu_count() or 1
    max_workers = max(1, min(int(max_workers), len(pairs)))
    chunk_size = cfg.get("filter_exists_chunk_size")
    if chunk_size is None:
        chunk_size = max(1, len(pairs) // max_workers)
    chunk_size = max(1, int(chunk_size))

    if use_threaded and max_workers > 1 and len(pairs) > chunk_size:
        try:
            return _run_threaded_exists_check_break_both(
                pairs, sp, dist, w, n, t, max_workers, chunk_size
            )
        except Exception:
            pass  # Fall through to sequential

    return _exists_violation_sequential(pairs, sp, dist, w, n, t)


def _filter_fallback_reinsert_predicate(
    selected: list[tuple[int, int]],
    without: list[tuple[int, int]],
    i: int,
    j: int,
    dist: np.ndarray,
    w: np.ndarray,
    n: int,
    t: float,
    config: Optional[dict[str, Any]],
) -> bool:
    """
    SP-50/SP-57: Global exists (r,s) predicate. Called only when local check
    delta(pq)*w(pq) > t*|pq| is false. Returns True to reinsert iff any
    unordered pair (r,s) in full-point-set domain violates delta(rs)*w(rs) > t*|rs|
    in graph without. Fallback logic is independent from E_input subset semantics.
    Break-both: stop on first violation.
    """
    _ = selected, i, j
    return _exists_violation_pairs(without, dist, w, n, t, config)


def _filter_pass(
    selected: list[tuple[int, int]],
    dist: np.ndarray,
    w: np.ndarray,
    n: int,
    t: float,
    config: Optional[dict[str, Any]] = None,
) -> tuple[list[tuple[int, int]], int]:
    """
    One filter pass: remove/reinsert validation.
    Strict sequence: local-first then all pairs validation. Never call all pairs validation
    when local check is true.

    Decision sequence:
    1) Remove (p,q).
    2) Evaluate local check: delta(pq)*w(pq) > t*|pq| (on graph without (p,q)).
    3) If local check true → reinsert (fallback NOT called).
    4) If local check false → call fallback predicate; reinsert if it returns true.
    5) Otherwise keep removed.

    Iteration order: descending by length.

    Uses incremental remove-edge path when apsp_mode != "full"; full APSP
    fallback for A/B validation and rollback.
    """
    if len(selected) <= 1:
        return selected, 0

    apsp_mode = resolve_apsp_mode(config)
    use_incremental = apsp_mode != APSP_MODE_FULL

    # Build ordered edge list
    sel_arr = np.array(selected, dtype=np.int64)
    lengths = dist[sel_arr[:, 0], sel_arr[:, 1]].astype(np.float64)
    order = np.argsort(-lengths)  # descending
    ordered = [selected[i] for i in order]

    removed = 0
    for i, j in progress_iter(
        ordered,
        total=len(ordered),
        desc="filter:pass",
        config=config,
    ):
        edge_canon = tuple(sorted((i, j)))
        # 1) Remove (p,q)
        without = [(a, b) for (a, b) in selected if tuple(sorted((a, b))) != edge_canon]
        if len(without) == len(selected):
            continue  # edge not in selected (canonicalization mismatch)

        if use_incremental:
            E_curr = np.array(selected, dtype=np.int64)
            delta_E_prime = apsp_shortest_path_after_removal(E_curr, i, j, n, dist)
        else:
            E_without = np.array(without, dtype=np.int64)
            sp = _all_pairs_shortest_paths(E_without, n, dist)
            delta_E_prime = sp[i, j]

        pq = dist[i, j]
        wij = w[i, j]

        # 2) Local check: delta(pq)*w(pq) > t*|pq| → reinsert
        if np.isinf(delta_E_prime):
            # Disconnected: treat as violated, reinsert
            continue
        if wij * delta_E_prime > t * pq:
            # Local check true → reinsert (keep in selected)
            continue

        # 3) Local check false → call fallback predicate
        if _filter_fallback_reinsert_predicate(selected, without, i, j, dist, w, n, t, config):
            # Fallback says reinsert
            continue

        # 4) Keep removed
        selected = without
        removed += 1

    return selected, removed


@register("filter")
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
    Filter algorithm: remove/reinsert validation.


    (filter): Remove/reinsert validation. Iterate edges descending
        by length. Remove iff w(p,q)*delta_E' <= t*|pq| (graph valid without it).

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

    return np.array(_filter_pass(candidates.tolist(), dist, w, n, t, config=config)[0], dtype=np.int64)
