"""
Faster DGF (Descending Greedy Filter) v7
=========================================
Batched APSP with full recursive bisection. No precheck, no deferral.

Edges sorted by descending weight are buffered into batches. Each batch is
tested with one APSP + violation check. If clean, all removed in one shot.
If violated, recurse on both halves to find all necessary edges.

Interface:
    faster_dgf_one_pass(edges, dist_np, n, t) -> (remaining_edges, n_removed)
"""

from __future__ import annotations

import os
import time

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

# Pre-compiled FW kernel — one launch per k-iteration instead of 2 CuPy ops
_fw_kernel = None
if _CUDA:
    _fw_kernel = cp.RawKernel(r'''
    extern "C" __global__
    void fw_step(double* d, int n, int k) {
        int idx = blockDim.x * blockIdx.x + threadIdx.x;
        int i = idx / n;
        int j = idx % n;
        if (i < n && j < n) {
            double newval = d[i * n + k] + d[k * n + j];
            if (newval < d[i * n + j])
                d[i * n + j] = newval;
        }
    }
    ''', 'fw_step')

try:
    from native_dgf import dgf_one_pass_v2 as native_dgf_one_pass_v2
    _NATIVE_V2 = True
except Exception:
    _NATIVE_V2 = False

try:
    from native_dgf import dgf_one_pass_native
    _NATIVE_V1 = True
except Exception:
    _NATIVE_V1 = False

_DGF_NATIVE_ENABLED = os.environ.get("DGF_NATIVE", "1") != "0"


def _log(msg: str) -> None:
    tqdm.write(f"    [dgf] {msg}")


# ── CSR helpers ───────────────────────────────────────────────────────────────

def _build_csr(exists: np.ndarray, dist_np: np.ndarray, n: int) -> sp_sparse.csr_matrix:
    """Build sorted CSR from boolean exists matrix."""
    rows, cols = np.where(exists)
    data = dist_np[rows, cols]
    csr = sp_sparse.csr_matrix((data, (rows, cols)), shape=(n, n), dtype=np.float64)
    csr.sort_indices()
    return csr


def _csr_find(csr: sp_sparse.csr_matrix, row: int, col: int) -> int:
    """Find data index for (row,col). Returns -1 if absent."""
    start, end = csr.indptr[row], csr.indptr[row + 1]
    pos = np.searchsorted(csr.indices[start:end], col)
    if pos < end - start and csr.indices[start + pos] == col:
        return start + pos
    return -1


def _csr_remove(csr: sp_sparse.csr_matrix, u: int, v: int) -> None:
    """Set (u,v) and (v,u) to inf (Dijkstra ignores inf-weight edges). O(log degree)."""
    for a, b in [(u, v), (v, u)]:
        idx = _csr_find(csr, a, b)
        if idx >= 0:
            csr.data[idx] = np.inf


def _csr_restore(csr: sp_sparse.csr_matrix, u: int, v: int, w: float) -> None:
    """Restore (u,v) and (v,u). O(log degree)."""
    for a, b in [(u, v), (v, u)]:
        idx = _csr_find(csr, a, b)
        if idx >= 0:
            csr.data[idx] = w


# ── Violation check ──────────────────────────────────────────────────────────

def _any_violation(sp_matrix: np.ndarray, dist_np: np.ndarray, t: float, n: int) -> bool:
    """Check if any pair (i,j) has sp[i,j] > t * dist[i,j]. GPU-accelerated if available."""
    if _CUDA:
        sp_gpu = cp.asarray(sp_matrix)
        dist_gpu = cp.asarray(dist_np[:n, :n])
        violated = cp.any(sp_gpu > t * dist_gpu + 1e-9)
        return bool(violated)
    else:
        return bool(np.any(sp_matrix > t * dist_np[:n, :n] + 1e-9))


# ── Parallel APSP ────────────────────────────────────────────────────────────

def _get_n_workers() -> int:
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        return os.cpu_count() or 1


def _parallel_apsp(csr: sp_sparse.csr_matrix, n: int, stats: dict | None = None) -> np.ndarray:
    """Multi-threaded APSP. GPU FW for moderate n, else parallel CPU Dijkstra."""
    nw = _get_n_workers()

    # GPU Floyd-Warshall with RawKernel: 1 kernel launch per k (not 2 CuPy ops)
    if _CUDA and _fw_kernel is not None and n <= 6000:
        t0 = time.perf_counter()
        # Convert CSR to dense adjacency
        adj = np.full((n, n), np.inf, dtype=np.float64)
        np.fill_diagonal(adj, 0.0)
        rows, cols = csr.nonzero()
        vals = np.array(csr[rows, cols]).ravel()
        # Only copy finite values (inf entries stay as inf in adj)
        finite_mask = np.isfinite(vals)
        if np.any(finite_mask):
            adj[rows[finite_mask], cols[finite_mask]] = vals[finite_mask]
        # GPU FW via RawKernel
        d = cp.asarray(adj)
        threads = 256
        blocks = (n * n + threads - 1) // threads
        n_i32 = np.int32(n)
        for k in range(n):
            _fw_kernel((blocks,), (threads,), (d, n_i32, np.int32(k)))
        result = cp.asnumpy(d)
        elapsed = time.perf_counter() - t0
        if stats is not None:
            stats["apsp_calls"] = stats.get("apsp_calls", 0) + 1
            stats["apsp_total_s"] = stats.get("apsp_total_s", 0.0) + elapsed
            stats["apsp_last_s"] = elapsed
            stats["apsp_method"] = "gpu_fw_rawkernel"
            stats["apsp_workers"] = 1
        return result

    if n < 200 or nw <= 1:
        t0 = time.perf_counter()
        result = sp_shortest_path(csr, directed=False)
        elapsed = time.perf_counter() - t0
        if stats is not None:
            stats["apsp_calls"] = stats.get("apsp_calls", 0) + 1
            stats["apsp_total_s"] = stats.get("apsp_total_s", 0.0) + elapsed
            stats["apsp_last_s"] = elapsed
            stats["apsp_method"] = "sequential"
            stats["apsp_workers"] = 1
        return result

    chunk_size = min(64, (n + nw - 1) // nw)  # small chunks for dynamic load balancing
    chunks = [list(range(i, min(i + chunk_size, n))) for i in range(0, n, chunk_size)]

    chunk_times: list[float] = []

    def _run(indices):
        tc = time.perf_counter()
        result = sp_dijkstra(csr, indices=indices, directed=False)
        chunk_times.append(time.perf_counter() - tc)
        return result

    t0 = time.perf_counter()
    with ThreadPoolExecutor(max_workers=nw) as ex:
        results = list(ex.map(_run, chunks))
    elapsed = time.perf_counter() - t0

    if stats is not None:
        stats["apsp_calls"] = stats.get("apsp_calls", 0) + 1
        stats["apsp_total_s"] = stats.get("apsp_total_s", 0.0) + elapsed
        stats["apsp_last_s"] = elapsed
        stats["apsp_method"] = "parallel"
        stats["apsp_workers"] = nw
        stats["apsp_chunks"] = len(chunks)
        stats["apsp_chunk_sizes"] = [len(c) for c in chunks]
        if chunk_times:
            stats["apsp_chunk_min_s"] = min(chunk_times)
            stats["apsp_chunk_max_s"] = max(chunk_times)
            stats["apsp_chunk_avg_s"] = sum(chunk_times) / len(chunk_times)

    return np.vstack(results)


# ── Binary search ────────────────────────────────────────────────────────────

def _bisect(
    csr: sp_sparse.csr_matrix,
    exists: np.ndarray,
    dist_np: np.ndarray,
    n: int,
    t: float,
    edges: list[tuple[int, int, float]],
    apsp_stats: dict,
    bisect_stats: dict,
    edge_bar: tqdm,
    t_total: float,
    n_edges_total: int,
) -> int:
    """
    All edges in *edges* are already removed from csr.
    Sorted in descending weight order.

    - No violation → all removed (1→2 together).
    - Single edge with violation → necessary (1→3).
    - Otherwise: restore all, recurse on both halves (longer first).

    Returns number of edges removed.
    """
    if not edges:
        return 0

    sp_matrix = _parallel_apsp(csr, n, apsp_stats)

    if not _any_violation(sp_matrix, dist_np, t, n):
        # All removable
        for p, q, _w in edges:
            exists[p, q] = False
            exists[q, p] = False
        edge_bar.update(len(edges))
        bisect_stats["removed"] += len(edges)
        # Log when a large chunk is cleared
        if len(edges) >= 64:
            _log(f"cleared {len(edges)} edges | "
                 f"removed={bisect_stats['removed']} "
                 f"necessary={bisect_stats['necessary']} "
                 f"settled={bisect_stats['removed']+bisect_stats['necessary']}/{n_edges_total} | "
                 f"apsp: {apsp_stats.get('apsp_calls', 0)} calls, "
                 f"{apsp_stats.get('apsp_total_s', 0):.1f}s | "
                 f"elapsed={time.perf_counter()-t_total:.1f}s")
        return len(edges)

    # Single edge → necessary
    if len(edges) == 1:
        p, q, w = edges[0]
        _csr_restore(csr, p, q, w)
        bisect_stats["necessary"] += 1
        edge_bar.update(1)
        return 0

    # Restore all, split, recurse on both halves
    for p, q, w in edges:
        _csr_restore(csr, p, q, w)

    mid = len(edges) // 2
    first_half, second_half = edges[:mid], edges[mid:]

    # First half: longer/heavier edges
    for p, q, _w in first_half:
        _csr_remove(csr, p, q)
    first_removed = _bisect(csr, exists, dist_np, n, t, first_half,
                            apsp_stats, bisect_stats, edge_bar, t_total, n_edges_total)

    # Second half: shorter/lighter edges
    for p, q, _w in second_half:
        _csr_remove(csr, p, q)
    second_removed = _bisect(csr, exists, dist_np, n, t, second_half,
                             apsp_stats, bisect_stats, edge_bar, t_total, n_edges_total)

    return first_removed + second_removed


# ── Main ─────────────────────────────────────────────────────────────────────

def faster_dgf_one_pass(
    edges: list[tuple[int, int]],
    dist_np: np.ndarray,
    n: int,
    t: float,
) -> tuple[list[tuple[int, int]], int]:
    """
    Descending Greedy Filter v8 — global binary search over all edges.
    No batches, no precheck, no deferral.
    """
    if len(edges) <= 1:
        return list(edges), 0

    t_total = time.perf_counter()
    sel = sorted(edges, key=lambda e: -dist_np[e[0], e[1]])
    n_edges = len(sel)

    _log(f"n={n}, edges={n_edges}, t={t}, CUDA={_CUDA}, "
         f"native_v2={_NATIVE_V2}, native_v1={_NATIVE_V1}")

    # ── Native Rust fast path ─────────────────────────────────────────
    if _DGF_NATIVE_ENABLED:
        if _NATIVE_V2:
            _log("trying native Rust v2 (bidirectional A*)")
            try:
                t0 = time.perf_counter()
                result = native_dgf_one_pass_v2(sel, dist_np, n, t)
                _log(f"native v2 done: {time.perf_counter()-t0:.1f}s")
                return result
            except Exception as e:
                _log(f"native v2 FAILED: {e}, falling back")
        if _NATIVE_V1:
            _log("trying native Rust v1")
            try:
                t0 = time.perf_counter()
                result = dgf_one_pass_native(sel, dist_np, n, t)
                _log(f"native v1 done: {time.perf_counter()-t0:.1f}s")
                return result
            except Exception as e:
                _log(f"native v1 FAILED: {e}, falling back")

    # ── exists matrix + CSR ───────────────────────────────────────────
    exists = np.zeros((n, n), dtype=bool)
    for u, v in sel:
        exists[u, v] = True
        exists[v, u] = True

    csr = _build_csr(exists, dist_np, n)

    apsp_stats: dict = {}
    bisect_stats: dict = {"removed": 0, "necessary": 0}

    _log(f"global binary search, CUDA={'yes' if _CUDA else 'no'}")

    # ── Build edge list with weights, remove all from CSR ─────────────
    all_edges = [(p, q, dist_np[p, q]) for p, q in sel]
    for p, q, _w in all_edges:
        _csr_remove(csr, p, q)

    edge_bar = tqdm(total=n_edges, desc="dgf", unit="edge", leave=False)

    # ── One global binary search ─────────────────────────────────────
    removed = _bisect(csr, exists, dist_np, n, t, all_edges,
                      apsp_stats, bisect_stats, edge_bar, t_total, n_edges)

    edge_bar.close()

    # ── Final report ──────────────────────────────────────────────────
    total_elapsed = time.perf_counter() - t_total
    _log(f"DONE in {total_elapsed:.1f}s | "
         f"removed={removed} necessary={bisect_stats['necessary']} "
         f"settled={removed + bisect_stats['necessary']}/{n_edges} | "
         f"apsp: {apsp_stats.get('apsp_calls', 0)} calls, "
         f"{apsp_stats.get('apsp_total_s', 0):.2f}s")

    remaining = [e for e in sel if exists[e[0], e[1]]]
    return remaining, removed
