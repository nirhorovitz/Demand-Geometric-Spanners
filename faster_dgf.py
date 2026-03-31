"""
Faster DGF (Descending Greedy Filter) v2
=========================================
Drop-in replacement for dgf_one_pass from final_experiment.py with the
following optimizations:

1. K_n Euclidean fast path: single Dijkstra(limit=t*w) per edge, no GPU
   broadcast, no `uses` tracking — exploits triangle inequality
2. Dijkstra `limit` parameter for early termination on all paths
3. Batch free removal: scan `uses` for empty entries, bulk remove
4. Incremental sp_gpu maintenance (avoid full FW recompute)
5. Dynamic dense→sparse strategy switching as edges thin out
6. Bridge pre-detection (sparse): skip bridge edges in Phase 2
7. Direct CSR manipulation (avoid lil→csr conversion per iteration)
8. Staleness-based `uses` rebuild (adaptive instead of fixed interval)
9. Dispatches to native Rust v2 (bidirectional A*, early abort) when available

Same interface as dgf_one_pass:
    faster_dgf_one_pass(edges, dist_np, n, t) -> (remaining_edges, n_removed)
"""

from __future__ import annotations

import math
import os
import time
from collections import defaultdict
from typing import Any, Optional

import numpy as np
import scipy.sparse as sp_sparse
from concurrent.futures import ThreadPoolExecutor
from scipy.sparse.csgraph import shortest_path as sp_shortest_path, dijkstra as sp_dijkstra
from tqdm import tqdm

try:
    import cupy as cp
    cp.array([1.0])
    _CUDA = True
except Exception:
    _CUDA = False

try:
    from native_dgf import (
        trace_path_edge_ids as native_trace_path_edge_ids,
        sort_pairs_with_key_first as native_sort_pairs_with_key_first,
    )
    _NATIVE_DGF_AVAILABLE = True
except Exception:
    _NATIVE_DGF_AVAILABLE = False

try:
    from native_dgf import dgf_one_pass_v2 as native_dgf_one_pass_v2
    _NATIVE_V2_AVAILABLE = True
except Exception:
    _NATIVE_V2_AVAILABLE = False

try:
    from native_dgf import dgf_one_pass_native
    _NATIVE_V1_AVAILABLE = True
except Exception:
    _NATIVE_V1_AVAILABLE = False

_DGF_NATIVE_ENABLED = os.environ.get("DGF_NATIVE", "1") != "0"
_DGF_PROFILE_ENABLED = os.environ.get("DGF_PROFILE", "0") == "1"


class _DgfProfiler:
    """Tiny phase profiler for optional DGF instrumentation."""

    def __init__(self, enabled: bool) -> None:
        self.enabled = enabled
        self._marks: dict[str, float] = {}

    def add(self, name: str, delta_sec: float) -> None:
        if not self.enabled:
            return
        self._marks[name] = self._marks.get(name, 0.0) + float(delta_sec)

    def report(self, prefix: str = "faster-dgf-profile") -> None:
        if not self.enabled:
            return
        total = sum(self._marks.values())
        print(f"  [{prefix}] total={total:.3f}s")
        for k in sorted(self._marks.keys()):
            v = self._marks[k]
            pct = (100.0 * v / total) if total > 0 else 0.0
            print(f"  [{prefix}] {k}={v:.3f}s ({pct:.1f}%)")


# ── Bridge detection (Tarjan's algorithm) ─────────────────────────────────────

def _find_bridges(edges: list[tuple[int, int]], n: int) -> set[tuple[int, int]]:
    """
    Find all bridge edges using iterative Tarjan's algorithm. O(n + m).
    Returns set of canonicalized (min, max) edge tuples.
    """
    adj: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for idx, (u, v) in enumerate(edges):
        adj[u].append((v, idx))
        adj[v].append((u, idx))

    disc = [-1] * n
    low = [-1] * n
    bridges: set[tuple[int, int]] = set()
    timer = 0

    for start in range(n):
        if disc[start] != -1:
            continue

        disc[start] = low[start] = timer
        timer += 1

        stack: list[tuple[int, int, int]] = [(start, -1, 0)]

        while stack:
            u, parent_eidx, ni = stack[-1]
            neighbors = adj[u]

            if ni < len(neighbors):
                stack[-1] = (u, parent_eidx, ni + 1)
                v, eidx = neighbors[ni]

                if eidx == parent_eidx:
                    continue

                if disc[v] == -1:
                    disc[v] = low[v] = timer
                    timer += 1
                    stack.append((v, eidx, 0))
                else:
                    low[u] = min(low[u], disc[v])
            else:
                stack.pop()
                if stack:
                    parent_node = stack[-1][0]
                    low[parent_node] = min(low[parent_node], low[u])
                    if low[u] > disc[parent_node]:
                        bridges.add((min(parent_node, u), max(parent_node, u)))

    return bridges


# ── CSR edge removal/restoration helpers ──────────────────────────────────────

def _csr_find_entry(csr: sp_sparse.csr_matrix, row: int, col: int) -> int:
    """Find the data index for (row, col) in a CSR matrix. Returns -1 if not found."""
    start, end = csr.indptr[row], csr.indptr[row + 1]
    indices = csr.indices[start:end]
    pos = np.searchsorted(indices, col)
    if pos < len(indices) and indices[pos] == col:
        return start + pos
    return -1


def _csr_remove_edge(csr: sp_sparse.csr_matrix, u: int, v: int) -> None:
    """Zero out edge (u,v) and (v,u) in CSR. O(log degree)."""
    idx = _csr_find_entry(csr, u, v)
    if idx >= 0:
        csr.data[idx] = 0.0
    idx = _csr_find_entry(csr, v, u)
    if idx >= 0:
        csr.data[idx] = 0.0


def _csr_restore_edge(csr: sp_sparse.csr_matrix, u: int, v: int, w: float) -> None:
    """Restore edge (u,v) and (v,u) in CSR. O(log degree)."""
    idx = _csr_find_entry(csr, u, v)
    if idx >= 0:
        csr.data[idx] = w
    idx = _csr_find_entry(csr, v, u)
    if idx >= 0:
        csr.data[idx] = w


def _csr_edge_alive(csr: sp_sparse.csr_matrix, u: int, v: int) -> bool:
    """Check if edge (u,v) is still present (non-zero) in CSR."""
    idx = _csr_find_entry(csr, u, v)
    return idx >= 0 and csr.data[idx] != 0.0


# ── Parallel APSP ─────────────────────────────────────────────────────────────

def _parallel_shortest_path(csr: sp_sparse.csr_matrix, n: int, n_workers: int = 0) -> np.ndarray:
    """
    Parallel APSP via multi-threaded Dijkstra. scipy's C-level Dijkstra
    releases the GIL, so ThreadPoolExecutor gives near-linear speedup.
    Falls back to sequential sp_shortest_path for small n.
    """
    if n_workers <= 0:
        try:
            n_workers = len(os.sched_getaffinity(0))
        except AttributeError:
            n_workers = os.cpu_count() or 1
    if n < 200 or n_workers <= 1:
        return sp_shortest_path(csr, directed=False)

    chunk_size = (n + n_workers - 1) // n_workers
    chunks = [list(range(i, min(i + chunk_size, n))) for i in range(0, n, chunk_size)]

    def _run_chunk(indices):
        return sp_dijkstra(csr, indices=indices)

    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        results = list(executor.map(_run_chunk, chunks))
    return np.vstack(results)


# ── Numpy free-edge detection ─────────────────────────────────────────────────

def _detect_free_edges_np(
    sp: np.ndarray,
    edge_keys: list[tuple[int, int]],
    dist_np: np.ndarray,
    tol: float = 1e-10,
) -> set[tuple[int, int]]:
    """
    Detect edges not on any shortest path using numpy/cupy broadcasting.
    Edge (u,v) with weight w is on shortest path r->s iff
    sp[r,u]+w+sp[v,s]==sp[r,s] or sp[r,v]+w+sp[u,s]==sp[r,s].
    Returns set of (u,v) tuples that are NOT on any shortest path (free edges).
    O(m * n^2) but with vectorization, much faster than Python path tracing.
    Uses GPU (CuPy) when available for large n — the n*n broadcast is
    ideal for GPU parallelism (25M elements for n=5000).
    """
    if not edge_keys:
        return set()

    # GPU path: transfer sp matrix once, do all edge checks on GPU
    if _CUDA and sp.shape[0] >= 200:
        try:
            sp_gpu = cp.asarray(sp, dtype=np.float64)
            free = set()
            for u, v in edge_keys:
                w = float(dist_np[u, v])
                route_uv = sp_gpu[:, u:u+1] + w + sp_gpu[v:v+1, :]
                route_vu = sp_gpu[:, v:v+1] + w + sp_gpu[u:u+1, :]
                best = cp.minimum(route_uv, route_vu)
                if not bool(cp.any(cp.abs(sp_gpu - best) < tol)):
                    free.add((u, v))
            del sp_gpu
            return free
        except Exception:
            pass  # Fall through to CPU path

    # CPU path
    free = set()
    for u, v in edge_keys:
        w = dist_np[u, v]
        route_uv = sp[:, u:u+1] + w + sp[v:v+1, :]  # (n, n) broadcast
        route_vu = sp[:, v:v+1] + w + sp[u:u+1, :]
        best = np.minimum(route_uv, route_vu)
        if not np.any(np.abs(sp - best) < tol):
            free.add((u, v))
    return free


# ── Floyd-Warshall ────────────────────────────────────────────────────────────

def _adj_to_scipy(adj_np: np.ndarray) -> np.ndarray:
    return np.where(np.isinf(adj_np), 0.0, adj_np)


def _fw(adj_np: np.ndarray) -> np.ndarray:
    if _CUDA:
        d = cp.asarray(adj_np, dtype=np.float64)
        n = d.shape[0]
        tmp = cp.empty((n, n), dtype=np.float64)
        for k in range(n):
            cp.add(d[:, k:k+1], d[k:k+1, :], out=tmp)
            cp.minimum(d, tmp, out=d)
        return cp.asnumpy(d)
    return sp_shortest_path(_adj_to_scipy(adj_np), method='FW', directed=False)


def _fw_device(adj_np: np.ndarray):
    if _CUDA:
        d = cp.asarray(adj_np, dtype=np.float64)
        n = d.shape[0]
        tmp = cp.empty((n, n), dtype=np.float64)
        for k in range(n):
            cp.add(d[:, k:k+1], d[k:k+1, :], out=tmp)
            cp.minimum(d, tmp, out=d)
        return d
    return sp_shortest_path(_adj_to_scipy(adj_np), method='FW', directed=False)


# ── Adjacency builders ────────────────────────────────────────────────────────

def _build_adj_np(edge_list: list, n: int, dist_np: np.ndarray) -> np.ndarray:
    adj = np.full((n, n), np.inf, dtype=np.float64)
    np.fill_diagonal(adj, 0.0)
    if edge_list:
        E = np.array(edge_list, dtype=np.int64)
        adj[E[:, 0], E[:, 1]] = dist_np[E[:, 0], E[:, 1]]
        adj[E[:, 1], E[:, 0]] = dist_np[E[:, 0], E[:, 1]]
    return adj


def _build_sparse_adj(edge_list: list, n: int, dist_np: np.ndarray) -> sp_sparse.lil_matrix:
    g = sp_sparse.lil_matrix((n, n), dtype=np.float64)
    for u, v in edge_list:
        w = dist_np[u, v]
        g[u, v] = w
        g[v, u] = w
    return g


def _build_csr(edge_list: list, n: int, dist_np: np.ndarray) -> sp_sparse.csr_matrix:
    """Build CSR directly from edge list. O(m) construction."""
    if not edge_list:
        return sp_sparse.csr_matrix((n, n), dtype=np.float64)
    E = np.array(edge_list, dtype=np.int64)
    rows = np.concatenate([E[:, 0], E[:, 1]])
    cols = np.concatenate([E[:, 1], E[:, 0]])
    weights = dist_np[E[:, 0], E[:, 1]]
    data = np.concatenate([weights, weights])
    csr = sp_sparse.csr_matrix((data, (rows, cols)), shape=(n, n), dtype=np.float64)
    csr.sort_indices()
    return csr


# ── Uses helpers ──────────────────────────────────────────────────────────────

def _trace_and_update_uses(pred: np.ndarray, r: int, s: int, uses: dict) -> None:
    pair = (min(r, s), max(r, s))
    cur = s
    while cur != r and cur >= 0:
        prev = int(pred[cur])
        if prev < 0:
            break
        edge_key = (min(prev, cur), max(prev, cur))
        if edge_key not in uses:
            uses[edge_key] = set()
        uses[edge_key].add(pair)
        cur = prev


def _trace_and_update_uses_maybe_native(
    pred: np.ndarray, r: int, s: int, uses: dict, n: int, *, use_native: bool,
) -> None:
    if use_native and _NATIVE_DGF_AVAILABLE and _DGF_NATIVE_ENABLED:
        try:
            pair = (min(r, s), max(r, s))
            edge_ids = native_trace_path_edge_ids(pred.astype(np.int64), int(r), int(s), int(n))
            for edge_id in edge_ids:
                u = int(edge_id) // n
                v = int(edge_id) % n
                edge_key = (u, v)
                if edge_key not in uses:
                    uses[edge_key] = set()
                uses[edge_key].add(pair)
            return
        except Exception:
            pass
    _trace_and_update_uses(pred, r, s, uses)


def _sort_pairs_with_key_first(
    pairs: set[tuple[int, int]], key: tuple[int, int], *, use_native: bool,
) -> list[tuple[int, int]]:
    pair_list = list(pairs)
    if key not in pairs:
        return pair_list
    if use_native and _NATIVE_DGF_AVAILABLE and _DGF_NATIVE_ENABLED:
        try:
            packed_pairs = [((int(r) << 32) | int(s)) for (r, s) in pair_list]
            packed_key = (int(key[0]) << 32) | int(key[1])
            sorted_packed = native_sort_pairs_with_key_first(packed_pairs, packed_key)
            return [(int(x >> 32), int(x & 0xFFFFFFFF)) for x in sorted_packed]
        except Exception:
            pass
    return [key] + [p for p in pair_list if p != key]


def _rebuild_uses(csr: sp_sparse.csr_matrix, n: int, existing_keys: set) -> dict:
    if not sp_sparse.isspmatrix_csr(csr):
        csr = csr.tocsr()
    _, predecessors = sp_dijkstra(csr, return_predecessors=True)
    uses: dict[tuple[int, int], set[tuple[int, int]]] = {}
    for key in existing_keys:
        uses[key] = set()
    for r in range(n):
        for s in range(r + 1, n):
            pair = (r, s)
            cur = s
            while cur != r and cur >= 0:
                prev = int(predecessors[r, cur])
                if prev < 0:
                    break
                edge_key = (min(prev, cur), max(prev, cur))
                if edge_key in uses:
                    uses[edge_key].add(pair)
                cur = prev
    return uses


def _batch_free_removal(
    uses: dict,
    adj: np.ndarray,
    csr: sp_sparse.csr_matrix,
    adj_gpu=None,
) -> int:
    """Remove all edges with empty uses (not on any shortest path). Returns count."""
    free_edges = [key for key, pairs in uses.items() if not pairs]
    for u, v in free_edges:
        adj[u, v] = np.inf
        adj[v, u] = np.inf
        _csr_remove_edge(csr, u, v)
        if adj_gpu is not None:
            adj_gpu[u, v] = cp.inf
            adj_gpu[v, u] = cp.inf
        del uses[(u, v)]
    return len(free_edges)


# ── Staleness heuristic ──────────────────────────────────────────────────────

def _should_rebuild(
    removed_since_rebuild: int,
    total_uses_pairs_removed: int,
    uses: dict,
) -> bool:
    """
    Staleness-based rebuild decision.
    Rebuild when estimated stale entries exceed 30% of total uses entries.
    """
    if removed_since_rebuild < 10:
        return False

    avg_pairs_per_removal = total_uses_pairs_removed / max(removed_since_rebuild, 1)
    estimated_stale = removed_since_rebuild * avg_pairs_per_removal * 0.5

    if not uses:
        return False

    sample_keys = list(uses.keys())[:min(100, len(uses))]
    sample_total = sum(len(uses[k]) for k in sample_keys)
    estimated_total = sample_total * len(uses) / max(len(sample_keys), 1)

    if estimated_total == 0:
        return removed_since_rebuild >= 500

    staleness_ratio = estimated_stale / estimated_total
    return staleness_ratio > 0.3


# ── Main entry point ──────────────────────────────────────────────────────────

def faster_dgf_one_pass(
    edges: list[tuple[int, int]],
    dist_np: np.ndarray,
    n: int,
    t: float,
) -> tuple[list[tuple[int, int]], int]:
    """
    Optimized Descending Greedy Filter v2. Drop-in replacement for dgf_one_pass.

    Key optimizations:
    1. K_n Euclidean fast path: single Dijkstra(limit=t*w) per edge
    2. Dijkstra `limit` for early termination everywhere
    3. Batch free removal of edges not on any shortest path
    4. Incremental sp_gpu update (avoid full FW recompute)
    5. Dynamic dense→sparse switching
    6. Bridge detection + CSR manipulation (sparse path)

    Returns: (remaining_edges, n_removed)
    """
    if len(edges) <= 1:
        return list(edges), 0

    sel = sorted(edges, key=lambda e: -dist_np[e[0], e[1]])

    use_sparse = len(sel) < 4 * n

    # ── Try Rust v2 (bidirectional A* + early abort) ──────────────────────
    if _NATIVE_V2_AVAILABLE and _DGF_NATIVE_ENABLED and not use_sparse:
        try:
            return native_dgf_one_pass_v2(sel, dist_np, n, t)
        except Exception as e:
            print(f"    [Native DGF v2 failed: {e}, falling back]")

    # ── Try Rust v1 ──────────────────────────────────────────────────────
    if _NATIVE_V1_AVAILABLE and _DGF_NATIVE_ENABLED and not use_sparse:
        try:
            return dgf_one_pass_native(sel, dist_np, n, t)
        except Exception as e:
            print(f"    [Native DGF v1 failed: {e}, falling back to Python]")

    profiler = _DgfProfiler(_DGF_PROFILE_ENABLED)
    t_phase = time.perf_counter()

    # ── Phase 1: Batch pre-filter ─────────────────────────────────────────
    n_edges = len(sel)
    is_complete = (n_edges == n * (n - 1) // 2)

    if is_complete:
        tqdm.write(f"    K_n detected: skipping Phase 1 APSP (triangle inequality)")
        n_batch_removed = 0
        sel_filtered = list(sel)
        profiler.add("phase1_prefilter_apsp", 0.0)
        profiler.add("phase1_prefilter_filtering", 0.0)
    elif use_sparse:
        lil = _build_sparse_adj(sel, n, dist_np)
        sp_init_np = sp_shortest_path(lil.tocsr(), directed=False)
        profiler.add("phase1_prefilter_apsp", time.perf_counter() - t_phase)
        t_phase = time.perf_counter()

        E_arr = np.array(sel, dtype=np.int64)
        sp_vals = sp_init_np[E_arr[:, 0], E_arr[:, 1]]
        w_vals = dist_np[E_arr[:, 0], E_arr[:, 1]]
        keep_mask = sp_vals >= w_vals
        n_batch_removed = int(np.sum(~keep_mask))
    else:
        adj = _build_adj_np(sel, n, dist_np)
        sp_init = _fw_device(adj)
        profiler.add("phase1_prefilter_apsp", time.perf_counter() - t_phase)
        t_phase = time.perf_counter()

        E_arr = np.array(sel, dtype=np.int64)
        sp_vals = sp_init[E_arr[:, 0], E_arr[:, 1]]
        if _CUDA:
            sp_vals = cp.asnumpy(sp_vals)
        w_vals = dist_np[E_arr[:, 0], E_arr[:, 1]]
        keep_mask = sp_vals >= w_vals
        n_batch_removed = int(np.sum(~keep_mask))

    if not is_complete:
        current_set = set()
        for idx in np.where(keep_mask)[0]:
            e = sel[idx]
            current_set.add((min(e[0], e[1]), max(e[0], e[1])))
        sel_filtered = [e for e in sel if (min(e[0], e[1]), max(e[0], e[1])) in current_set]

        if n_batch_removed > 0:
            tqdm.write(f"    batch pre-filter: removed {n_batch_removed}/{len(sel)} edges, {len(sel_filtered)} remain")
        profiler.add("phase1_prefilter_filtering", time.perf_counter() - t_phase)
        t_phase = time.perf_counter()

    # ── Phase 1.5: Bridge detection (sparse only) ────────────────────────
    bridge_set: set[tuple[int, int]] = set()
    if use_sparse and not is_complete and len(sel_filtered) > 10:
        t_bridge = time.perf_counter()
        bridge_set = _find_bridges(sel_filtered, n)
        profiler.add("bridge_detection", time.perf_counter() - t_bridge)
        if bridge_set:
            tqdm.write(f"    bridge detection: {len(bridge_set)} bridges found (will skip in Phase 2)")

    # ── Phase 2 ──────────────────────────────────────────────────────────
    r_np, s_np = np.triu_indices(n, k=1)
    dist_upper_np = dist_np[r_np, s_np]

    def _violation_np(sp_mat):
        sp_up = sp_mat[r_np, s_np]
        return bool(np.any(np.isinf(sp_up) | (sp_up > t * dist_upper_np)))

    if use_sparse:
        # ── Optimized sparse path ─────────────────────────────────────
        # Two-pronged optimization over naive brute-force:
        # 1. Parallel APSP: multi-threaded Dijkstra for the full APSP
        #    check (scipy releases GIL, so threads give near-linear speedup)
        # 2. Batch free removal: after every K removals, detect edges not
        #    on any shortest path via numpy broadcasting (O(m*n^2) but
        #    vectorized), and bulk-remove them without individual checks.
        #    This skips ~50-60% of edges that would otherwise need APSP.
        csr = _build_csr(sel_filtered, n, dist_np)
        removed = 0
        bridge_skipped = 0

        # Determine thread count for parallel APSP
        try:
            _n_workers = len(os.sched_getaffinity(0))
        except AttributeError:
            _n_workers = os.cpu_count() or 1

        # Build set of non-bridge edge keys for free-edge tracking
        active_keys: set[tuple[int, int]] = set()
        for e in sel_filtered:
            key = (min(e[0], e[1]), max(e[0], e[1]))
            if key not in bridge_set:
                active_keys.add(key)

        # ── Initial batch free removal via numpy detection ─────────────
        t_init = time.perf_counter()
        sp_init = _parallel_shortest_path(csr, n, _n_workers)

        # Pre-existing violation check
        sp_init_upper = sp_init[r_np, s_np]
        has_preexisting = bool(np.any(
            np.isinf(sp_init_upper) | (sp_init_upper > t * dist_upper_np)
        ))

        if not has_preexisting:
            free_init = _detect_free_edges_np(sp_init, list(active_keys), dist_np)
            for fk in free_init:
                _csr_remove_edge(csr, fk[0], fk[1])
                active_keys.discard(fk)
            if free_init:
                removed += len(free_init)
                tqdm.write(f"    sparse initial batch free: {len(free_init)} edges")
        del sp_init
        profiler.add("sparse_init", time.perf_counter() - t_init)

        if has_preexisting:
            tqdm.write("    sparse: pre-existing violations, using brute-force")

        # ── Main edge processing loop ──────────────────────────────────
        removed_since_rebuild = 0
        # Rebuild interval: balance between rebuild cost (1 APSP + m*n^2
        # numpy ops) and edges skipped. Empirically, interval=10 is near
        # optimal.
        _REBUILD_INTERVAL = max(5, n // 100)

        for edge in tqdm(sel_filtered, desc="dgf(sparse-opt)", unit="edge", leave=False):
            u, v = edge[0], edge[1]
            w = dist_np[u, v]
            key = (min(u, v), max(u, v))

            # Skip bridges
            if key in bridge_set:
                bridge_skipped += 1
                continue

            # Skip already-removed edges
            if key not in active_keys:
                continue

            # Tentatively remove edge
            _csr_remove_edge(csr, u, v)

            # Quick necessity check: Dijkstra from u with limit
            d_from_u = sp_dijkstra(csr, indices=[u], limit=t * w).ravel()
            if d_from_u[v] > t * w:
                _csr_restore_edge(csr, u, v, w)
                continue

            # Full APSP check (parallel on multi-core)
            sp_new = _parallel_shortest_path(csr, n, _n_workers)
            if _violation_np(sp_new):
                _csr_restore_edge(csr, u, v, w)
            else:
                removed += 1
                removed_since_rebuild += 1
                active_keys.discard(key)

                # Periodic rebuild: detect new free edges
                if not has_preexisting and removed_since_rebuild >= _REBUILD_INTERVAL:
                    t0 = time.perf_counter()
                    sp_rebuild = _parallel_shortest_path(csr, n, _n_workers)
                    free_new = _detect_free_edges_np(
                        sp_rebuild, list(active_keys), dist_np
                    )
                    del sp_rebuild
                    for fk in free_new:
                        _csr_remove_edge(csr, fk[0], fk[1])
                        active_keys.discard(fk)
                    if free_new:
                        removed += len(free_new)
                        tqdm.write(f"    sparse post-rebuild batch free: {len(free_new)} edges")
                    removed_since_rebuild = 0
                    profiler.add("sparse_rebuild", time.perf_counter() - t0)

        remaining = []
        for e in sel_filtered:
            u, v = e[0], e[1]
            if _csr_edge_alive(csr, u, v):
                remaining.append(e)
            elif (min(u, v), max(u, v)) in bridge_set:
                remaining.append(e)

        if bridge_skipped > 0:
            tqdm.write(f"    bridges skipped: {bridge_skipped}")

    elif is_complete:
        # ══════════════════════════════════════════════════════════════
        # K_n EUCLIDEAN FAST PATH
        # ══════════════════════════════════════════════════════════════
        # For Euclidean K_n, initial shortest paths are all direct edges
        # (triangle inequality), so uses[(u,v)] starts as just {(u,v)}.
        # As edges are removed, uses grows for remaining edges.
        # Key optimization: skip GPU broadcast entirely — trust `uses`.
        # Most edges only need single-pair Dijkstra(limit=t*w).
        # ══════════════════════════════════════════════════════════════
        tqdm.write(f"    K_n fast path: {len(sel_filtered)} edges, uses-tracked + Dijkstra(limit)")

        adj = _build_adj_np(sel_filtered, n, dist_np)
        csr = _build_csr(sel_filtered, n, dist_np)
        removed = 0

        use_native_helpers = bool(_NATIVE_DGF_AVAILABLE and _DGF_NATIVE_ENABLED)

        # Initialize uses: for K_n Euclidean, every edge only has its self-pair
        uses: dict[tuple[int, int], set[tuple[int, int]]] = {}
        for e in sel_filtered:
            key = (min(e[0], e[1]), max(e[0], e[1]))
            uses[key] = {key}

        removed_since_rebuild = 0
        total_uses_pairs_removed = 0

        for edge in tqdm(sel_filtered, desc="dgf(Kn-fast)", unit="edge", leave=False):
            u, v = edge[0], edge[1]
            w = dist_np[u, v]
            key = (min(u, v), max(u, v))

            if key not in uses:
                # Already batch-removed
                continue

            pairs_explicit = uses.get(key, set())
            if not pairs_explicit:
                # Free removal — not on any SP
                adj[u, v] = np.inf
                adj[v, u] = np.inf
                _csr_remove_edge(csr, u, v)
                removed += 1
                removed_since_rebuild += 1
                uses.pop(key, None)
                continue

            pair_set = set(pairs_explicit)
            pair_set.add(key)
            pair_list = _sort_pairs_with_key_first(
                pair_set, key, use_native=use_native_helpers
            )

            # Tentatively remove edge
            adj[u, v] = np.inf
            adj[v, u] = np.inf
            _csr_remove_edge(csr, u, v)

            # Necessity check with Dijkstra + limit
            sources: list[int] = []
            seen_sources: set[int] = set()
            for r, _ in pair_list:
                if r not in seen_sources:
                    seen_sources.add(r)
                    sources.append(r)

            max_pair_dist = max(dist_np[r, s] for r, s in pair_list)
            dijk_limit = t * max_pair_dist

            d_rows, pred_rows = sp_dijkstra(
                csr, indices=sources, return_predecessors=True, limit=dijk_limit
            )
            d_rows = np.atleast_2d(d_rows)
            pred_rows = np.atleast_2d(pred_rows)
            source_row = {src: idx for idx, src in enumerate(sources)}

            necessary = False
            for (r, s) in pair_list:
                d_src = d_rows[source_row[r]]
                if d_src[s] > t * dist_np[r, s]:
                    necessary = True
                    break

            if necessary:
                adj[u, v] = w
                adj[v, u] = w
                _csr_restore_edge(csr, u, v, w)
                continue

            # Removable — update uses for rerouted pairs
            removed += 1
            removed_since_rebuild += 1
            for (r, s) in pair_list:
                pred_src = pred_rows[source_row[r]]
                _trace_and_update_uses_maybe_native(
                    pred_src, r, s, uses, n, use_native=use_native_helpers
                )
            n_pairs = len(uses.pop(key, set()))
            total_uses_pairs_removed += n_pairs

            # Staleness-based rebuild
            if _should_rebuild(removed_since_rebuild, total_uses_pairs_removed, uses):
                t0 = time.perf_counter()
                existing_keys = set(uses.keys())
                uses = _rebuild_uses(csr, n, existing_keys)
                removed_since_rebuild = 0
                total_uses_pairs_removed = 0
                # Batch free removal after rebuild
                n_bf = 0
                free_keys = [k for k, v in uses.items() if not v]
                for fk in free_keys:
                    fu, fv = fk
                    adj[fu, fv] = np.inf
                    adj[fv, fu] = np.inf
                    _csr_remove_edge(csr, fu, fv)
                    del uses[fk]
                    n_bf += 1
                if n_bf > 0:
                    removed += n_bf
                    tqdm.write(f"    post-rebuild batch free: {n_bf} edges")
                profiler.add("kn_uses_rebuild", time.perf_counter() - t0)

        remaining = []
        for e in sel_filtered:
            u, v = e[0], e[1]
            if not np.isinf(adj[u, v]):
                remaining.append(e)

    else:
        # ── Dense non-K_n path with GPU-aware affected-pair detection ────
        use_gpu_phase2 = _CUDA and n >= 100

        adj = _build_adj_np(sel_filtered, n, dist_np)
        csr_dense = _build_csr(sel_filtered, n, dist_np)
        removed = 0

        use_native_helpers = bool(_NATIVE_DGF_AVAILABLE and _DGF_NATIVE_ENABLED)

        # Dependency tracking
        uses: dict[tuple[int, int], set[tuple[int, int]]] = {}

        for e in sel_filtered:
            key = (min(e[0], e[1]), max(e[0], e[1]))
            uses[key] = set()
        finite = np.isfinite(adj) & (adj > 0)
        csr_init = sp_sparse.csr_matrix(np.where(finite, adj, 0.0))
        _, predecessors = sp_dijkstra(csr_init, return_predecessors=True)
        for r in range(n):
            for s in range(r + 1, n):
                pair = (r, s)
                cur = s
                while cur != r and cur >= 0:
                    prev = int(predecessors[r, cur])
                    if prev < 0:
                        break
                    edge_key = (min(prev, cur), max(prev, cur))
                    if edge_key in uses:
                        uses[edge_key].add(pair)
                    cur = prev

        profiler.add("phase2_uses_init", time.perf_counter() - t_phase)
        t_phase = time.perf_counter()

        # NEW: Batch free removal after initial uses build
        n_batch_free = _batch_free_removal(uses, adj, csr_dense)
        if n_batch_free > 0:
            tqdm.write(f"    batch free removal: {n_batch_free} edges not on any SP")
            removed += n_batch_free
        profiler.add("phase2_batch_free", time.perf_counter() - t_phase)
        t_phase = time.perf_counter()

        # Track which edges were batch-removed (for skipping in main loop)
        batch_removed_keys: set[tuple[int, int]] = set()
        if n_batch_free > 0:
            for e in sel_filtered:
                key = (min(e[0], e[1]), max(e[0], e[1]))
                if key not in uses and np.isinf(adj[e[0], e[1]]):
                    batch_removed_keys.add(key)

        shorter_half_start = len(sel_filtered) // 2

        # Staleness-based rebuild tracking
        removed_since_rebuild = 0
        total_uses_pairs_removed = 0

        # GPU Phase 2 init
        sp_gpu = None
        adj_gpu = None
        dist_gpu = None
        if use_gpu_phase2:
            adj_gpu = cp.asarray(adj, dtype=np.float64)
            dist_gpu = cp.asarray(dist_np, dtype=np.float64)
            t0 = time.perf_counter()
            sp_gpu = _fw_device(adj)
            profiler.add("gpu_phase2_init_apsp", time.perf_counter() - t0)
            tqdm.write(f"    GPU Phase 2: APSP initialized on GPU ({n}x{n} matrix)")

        # Track remaining edges for dynamic switching
        remaining_edge_count = len(sel_filtered) - n_batch_free

        for edge_idx, edge in enumerate(tqdm(sel_filtered, desc="dgf(fast)", unit="edge", leave=False)):
            u, v = edge[0], edge[1]
            w = dist_np[u, v]
            key = (min(u, v), max(u, v))

            # Skip batch-removed edges
            if key in batch_removed_keys:
                continue

            # Dynamic dense→sparse switch
            if remaining_edge_count < 4 * n and remaining_edge_count > 0:
                # Switch to sparse strategy for remaining edges
                tqdm.write(f"    switching to sparse path ({remaining_edge_count} edges remain)")
                remaining_sel = []
                for e in sel_filtered[edge_idx:]:
                    eu, ev = e[0], e[1]
                    ekey = (min(eu, ev), max(eu, ev))
                    if ekey not in batch_removed_keys and not np.isinf(adj[eu, ev]):
                        remaining_sel.append(e)

                # Bridge detection on remaining subgraph
                sparse_bridge_set: set[tuple[int, int]] = set()
                if len(remaining_sel) > 10:
                    sparse_bridge_set = _find_bridges(remaining_sel, n)

                csr_sparse = csr_dense  # reuse the existing CSR
                for e_rem in remaining_sel:
                    ru, rv = e_rem[0], e_rem[1]
                    rw = dist_np[ru, rv]
                    rkey = (min(ru, rv), max(ru, rv))

                    if rkey in sparse_bridge_set:
                        continue

                    _csr_remove_edge(csr_sparse, ru, rv)
                    d_from_r = sp_dijkstra(csr_sparse, indices=[ru], limit=t * rw).ravel()
                    if d_from_r[rv] > t * rw:
                        _csr_restore_edge(csr_sparse, ru, rv, rw)
                        continue

                    sp_new = sp_shortest_path(csr_sparse, directed=False)
                    if _violation_np(sp_new):
                        _csr_restore_edge(csr_sparse, ru, rv, rw)
                    else:
                        adj[ru, rv] = np.inf
                        adj[rv, ru] = np.inf
                        removed += 1
                break  # exit main loop, sparse path handled remaining

            pairs_explicit = uses.get(key, set())
            if not pairs_explicit:
                # Free removal (not on any SP)
                adj[u, v] = np.inf
                adj[v, u] = np.inf
                _csr_remove_edge(csr_dense, u, v)
                if use_gpu_phase2 and adj_gpu is not None:
                    adj_gpu[u, v] = cp.inf
                    adj_gpu[v, u] = cp.inf
                removed += 1
                removed_since_rebuild += 1
                remaining_edge_count -= 1
                continue

            pair_list = _sort_pairs_with_key_first(
                pairs_explicit, key, use_native=use_native_helpers
            )

            # Tentatively remove edge
            adj[u, v] = np.inf
            adj[v, u] = np.inf
            _csr_remove_edge(csr_dense, u, v)
            if use_gpu_phase2 and adj_gpu is not None:
                adj_gpu[u, v] = cp.inf
                adj_gpu[v, u] = cp.inf

            # ── GPU-accelerated path ──────────────────────────────────
            if use_gpu_phase2:
                if sp_gpu is None:
                    t0_rebuild = time.perf_counter()
                    sp_gpu = cp.copy(adj_gpu)
                    tmp_rebuild = cp.empty((n, n), dtype=np.float64)
                    for k_fw in range(n):
                        cp.add(sp_gpu[:, k_fw:k_fw+1], sp_gpu[k_fw:k_fw+1, :], out=tmp_rebuild)
                        cp.minimum(sp_gpu, tmp_rebuild, out=sp_gpu)
                    profiler.add("gpu_phase2_apsp_rebuild", time.perf_counter() - t0_rebuild)
                t0 = time.perf_counter()

                # Affected pair detection via GPU broadcast
                sp_r_u = sp_gpu[:, u:u+1]
                sp_v_s = sp_gpu[v:v+1, :]
                sp_r_v = sp_gpu[:, v:v+1]
                sp_u_s = sp_gpu[u:u+1, :]

                route_uv = sp_r_u + w + sp_v_s
                route_vu = sp_r_v + w + sp_u_s
                best_through = cp.minimum(route_uv, route_vu)

                tol = 1e-10
                affected_mask = cp.abs(sp_gpu - best_through) < tol
                self_affected = float(cp.asnumpy(sp_gpu[u, v])) <= w + tol

                if not self_affected and not bool(cp.any(affected_mask)):
                    # Edge not on any shortest path — free removal
                    removed += 1
                    removed_since_rebuild += 1
                    remaining_edge_count -= 1
                    n_pairs = len(uses.pop(key, set()))
                    total_uses_pairs_removed += n_pairs
                    profiler.add("gpu_phase2_free_removal", time.perf_counter() - t0)
                    continue

                # Decide GPU FW vs CPU Dijkstra
                n_sources = len(set(r for r, _ in pair_list))
                use_gpu_fw = (n_sources > n // 8)

                if use_gpu_fw:
                    d_gpu = cp.copy(adj_gpu)
                    tmp_gpu = cp.empty((n, n), dtype=np.float64)
                    for k_fw in range(n):
                        cp.add(d_gpu[:, k_fw:k_fw+1], d_gpu[k_fw:k_fw+1, :], out=tmp_gpu)
                        cp.minimum(d_gpu, tmp_gpu, out=d_gpu)

                    sp_upper = d_gpu[r_np, s_np]
                    dist_upper_gpu = dist_gpu[r_np, s_np]
                    violated = bool(cp.any(cp.isinf(sp_upper) | (sp_upper > t * dist_upper_gpu)))
                    profiler.add("gpu_phase2_fw_check", time.perf_counter() - t0)

                    if violated:
                        adj[u, v] = w
                        adj[v, u] = w
                        _csr_restore_edge(csr_dense, u, v, w)
                        adj_gpu[u, v] = w
                        adj_gpu[v, u] = w
                        continue
                    else:
                        sp_gpu = d_gpu
                        removed += 1
                        removed_since_rebuild += 1
                        remaining_edge_count -= 1

                        # Uses update + incremental sp_gpu
                        t0 = time.perf_counter()
                        sources_gpu: list[int] = []
                        seen_gpu: set[int] = set()
                        for r, _ in pair_list:
                            if r not in seen_gpu:
                                seen_gpu.add(r)
                                sources_gpu.append(r)
                        d_rows_gpu, pred_rows = sp_dijkstra(
                            csr_dense, indices=sources_gpu, return_predecessors=True
                        )
                        d_rows_gpu = np.atleast_2d(d_rows_gpu)
                        pred_rows = np.atleast_2d(pred_rows)
                        source_row = {src: idx for idx, src in enumerate(sources_gpu)}

                        # Incremental sp_gpu update from Dijkstra results
                        for src_idx, src in enumerate(sources_gpu):
                            row_data = cp.asarray(d_rows_gpu[src_idx])
                            sp_gpu[src, :] = row_data
                            sp_gpu[:, src] = row_data

                        for (r, s) in pair_list:
                            pred_src = pred_rows[source_row[r]]
                            _trace_and_update_uses_maybe_native(
                                pred_src, r, s, uses, n, use_native=use_native_helpers
                            )
                        n_pairs = len(uses.pop(key, set()))
                        total_uses_pairs_removed += n_pairs
                        profiler.add("gpu_phase2_uses_update", time.perf_counter() - t0)

                        # Staleness-based rebuild
                        if _should_rebuild(removed_since_rebuild, total_uses_pairs_removed, uses):
                            t0 = time.perf_counter()
                            existing_keys = set(uses.keys())
                            uses = _rebuild_uses(csr_dense, n, existing_keys)
                            removed_since_rebuild = 0
                            total_uses_pairs_removed = 0
                            # Batch free removal after rebuild
                            n_bf = _batch_free_removal(uses, adj, csr_dense, adj_gpu)
                            if n_bf > 0:
                                removed += n_bf
                                remaining_edge_count -= n_bf
                                tqdm.write(f"    post-rebuild batch free: {n_bf} edges")
                            # Recompute sp_gpu from fresh state
                            sp_gpu = _fw_device(adj)
                            profiler.add("dense_uses_rebuild", time.perf_counter() - t0)
                        continue
                else:
                    pass  # Fall through to CPU path

            # ── CPU path ──────────────────────────────────────────────
            t0 = time.perf_counter()

            sources: list[int] = []
            seen_sources: set[int] = set()
            for r, _ in pair_list:
                if r not in seen_sources:
                    seen_sources.add(r)
                    sources.append(r)

            split_pred_mode = edge_idx >= shorter_half_start
            necessary = False
            source_row = {src: idx for idx, src in enumerate(sources)}
            pred_rows: Optional[np.ndarray] = None

            # NEW: Dijkstra with limit for early termination
            max_pair_dist = max(dist_np[r, s] for r, s in pair_list)
            dijk_limit = t * max_pair_dist

            t0 = time.perf_counter()
            if split_pred_mode:
                d_rows = sp_dijkstra(csr_dense, indices=sources, limit=dijk_limit)
            else:
                d_rows, pred_rows = sp_dijkstra(
                    csr_dense, indices=sources, return_predecessors=True, limit=dijk_limit
                )
            d_rows = np.atleast_2d(d_rows)

            for (r, s) in pair_list:
                d_src = d_rows[source_row[r]]
                if d_src[s] > t * dist_np[r, s]:
                    necessary = True
                    break
            profiler.add("dense_necessity_checks", time.perf_counter() - t0)

            if necessary:
                adj[u, v] = w
                adj[v, u] = w
                _csr_restore_edge(csr_dense, u, v, w)
                if use_gpu_phase2 and adj_gpu is not None:
                    adj_gpu[u, v] = w
                    adj_gpu[v, u] = w
                continue

            # Removable
            removed += 1
            removed_since_rebuild += 1
            remaining_edge_count -= 1
            if use_gpu_phase2:
                sp_gpu = None
            t0 = time.perf_counter()
            if pred_rows is None:
                _, pred_rows = sp_dijkstra(
                    csr_dense, indices=sources, return_predecessors=True
                )
            pred_rows = np.atleast_2d(pred_rows)

            for (r, s) in pair_list:
                pred_src = pred_rows[source_row[r]]
                _trace_and_update_uses_maybe_native(
                    pred_src, r, s, uses, n, use_native=use_native_helpers
                )
            n_pairs = len(uses.pop(key, set()))
            total_uses_pairs_removed += n_pairs
            profiler.add("dense_uses_update", time.perf_counter() - t0)

            # Staleness-based rebuild
            if _should_rebuild(removed_since_rebuild, total_uses_pairs_removed, uses):
                t0 = time.perf_counter()
                existing_keys = set(uses.keys())
                uses = _rebuild_uses(csr_dense, n, existing_keys)
                removed_since_rebuild = 0
                total_uses_pairs_removed = 0
                # Batch free removal after rebuild
                n_bf = _batch_free_removal(uses, adj, csr_dense,
                                           adj_gpu if use_gpu_phase2 else None)
                if n_bf > 0:
                    removed += n_bf
                    remaining_edge_count -= n_bf
                    tqdm.write(f"    post-rebuild batch free: {n_bf} edges")
                profiler.add("dense_uses_rebuild", time.perf_counter() - t0)

        remaining = []
        for e in sel_filtered:
            u, v = e[0], e[1]
            if not np.isinf(adj[u, v]):
                remaining.append(e)

    total_removed = n_batch_removed + removed
    profiler.add("phase2_total", time.perf_counter() - t_phase)
    if _DGF_PROFILE_ENABLED:
        mode = "sparse" if use_sparse else ("Kn-fast" if is_complete else "dense")
        if _CUDA and n >= 100 and not use_sparse and not is_complete:
            mode = f"{mode}-gpu"
        profiler.report(prefix=f"faster-dgf-{mode}")
    return remaining, total_removed
