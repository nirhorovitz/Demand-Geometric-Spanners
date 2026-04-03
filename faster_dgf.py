"""
Faster DGF (Descending Greedy Filter) v6
=========================================
Dijkstra precheck per edge, then batched APSP with binary-search bisection.

Edges that pass the cheap Dijkstra precheck are buffered. When the buffer
fills, one APSP + violation check tests the whole batch. If clean, all are
removed in one shot. If violated, recursive bisection finds which edges are
necessary — worst case 2x the per-edge cost, best case batch_size/1 speedup.

Interface:
    faster_dgf_one_pass(edges, dist_np, n, t) -> (remaining_edges, n_removed)
"""

from __future__ import annotations

import heapq
import os
import time
from typing import Optional

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


# ── Batch bisection ─────────────────────────────────────────────────────────

_BATCH_SIZE = int(os.environ.get("DGF_BATCH", "2056"))


def _bisect_batch(
    csr: sp_sparse.csr_matrix,
    exists: np.ndarray,
    dist_np: np.ndarray,
    n: int,
    t: float,
    batch: list[tuple[int, int, float]],
    apsp_stats: dict,
) -> tuple[int, list[tuple[int, int, float]]]:
    """
    All edges in *batch* are already removed from csr.
    Batch is in descending weight order: batch[:mid] = longer, batch[mid:] = shorter.

    On violation: restore the shorter half back to CSR and defer it (returned
    to caller so the main loop retries them later in normal order).  Recurse
    on the longer half only.

    Returns (n_removed, deferred_edges).
    """
    if not batch:
        return 0, []

    # One APSP for the whole batch
    sp_matrix = _parallel_apsp(csr, n, apsp_stats)
    if not _any_violation(sp_matrix, dist_np, t, n):
        # Entire batch removable
        for p, q, _w in batch:
            exists[p, q] = False
            exists[q, p] = False
        return len(batch), []

    # Single edge that causes violation → necessary, restore it
    if len(batch) == 1:
        p, q, w = batch[0]
        _csr_restore(csr, p, q, w)
        return 0, []

    # Split: longer = heavier edges (first half), shorter = lighter (second half)
    mid = len(batch) // 2
    longer, shorter = batch[:mid], batch[mid:]

    # Restore shorter half back to CSR — defer for later
    for p, q, w in shorter:
        _csr_restore(csr, p, q, w)

    # Recurse on longer half (still removed from CSR)
    longer_removed, longer_deferred = _bisect_batch(
        csr, exists, dist_np, n, t, longer, apsp_stats)

    return longer_removed, longer_deferred + shorter


# ── Main ─────────────────────────────────────────────────────────────────────

def faster_dgf_one_pass(
    edges: list[tuple[int, int]],
    dist_np: np.ndarray,
    n: int,
    t: float,
) -> tuple[list[tuple[int, int]], int]:
    """
    Descending Greedy Filter v6 — Dijkstra precheck + batched APSP bisection.
    """
    if len(edges) <= 1:
        return list(edges), 0

    t_total = time.perf_counter()
    sel = sorted(edges, key=lambda e: -dist_np[e[0], e[1]])
    n_edges = len(sel)

    _log(f"n={n}, edges={n_edges}, t={t}, CUDA={_CUDA}, "
         f"native_v2={_NATIVE_V2}, native_v1={_NATIVE_V1}, "
         f"batch_size={_BATCH_SIZE}")

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
    _log(f"precheck + batched APSP bisection, CUDA={'yes' if _CUDA else 'no'}")

    removed = 0
    precheck_necessary = 0
    precheck_passed = 0
    apsp_stats: dict = {}
    buffer: list[tuple[int, int, float]] = []

    # Max-heap: (-weight, p, q) so heaviest edges come first
    heap: list[tuple[float, int, int]] = [(-dist_np[e[0], e[1]], e[0], e[1]) for e in sel]
    heapq.heapify(heap)

    def _flush_buffer() -> int:
        nonlocal buffer
        if not buffer:
            return 0
        r, deferred = _bisect_batch(csr, exists, dist_np, n, t, buffer, apsp_stats)
        # Push shorter edges back onto heap — they'll be retried in order
        for dp, dq, dw in deferred:
            heapq.heappush(heap, (-dw, dp, dq))
        buffer = []
        return r

    # ── Precheck + batched APSP bisection ────────────────────────────
    processed = 0
    bar = tqdm(total=n_edges, desc="dgf", unit="edge", leave=False)

    while heap:
        neg_w, p, q = heapq.heappop(heap)
        w = -neg_w

        if not exists[p, q]:
            processed += 1
            bar.update(1)
            continue

        # 1. Tentatively remove
        _csr_remove(csr, p, q)

        # 2. Dijkstra precheck from p with limit = t*w
        d_pq = sp_dijkstra(csr, indices=[p], limit=t * w, directed=False).ravel()

        if d_pq[q] > t * w:
            # Necessary — restore immediately, never enters batch
            _csr_restore(csr, p, q, w)
            precheck_necessary += 1
            processed += 1
            bar.update(1)
            # Precheck failure signals we're in the critical zone —
            # flush whatever we've buffered so far
            if buffer:
                buf_size = len(buffer)
                r = _flush_buffer()
                removed += r
                _log(f"precheck-triggered flush: +{r}/{buf_size} removed "
                     f"(total {removed}) | "
                     f"apsp calls: {apsp_stats.get('apsp_calls', 0)} | "
                     f"elapsed: {time.perf_counter()-t_total:.1f}s")
            continue

        # 3. Passed precheck — edge stays removed in CSR, buffer it
        precheck_passed += 1
        buffer.append((p, q, w))
        processed += 1
        bar.update(1)

        if len(buffer) >= _BATCH_SIZE:
            r = _flush_buffer()
            removed += r
            _log(f"batch done: +{r} removed (total {removed}) | "
                 f"precheck: {precheck_necessary} necessary, "
                 f"{precheck_passed} passed | "
                 f"apsp calls: {apsp_stats.get('apsp_calls', 0)} | "
                 f"elapsed: {time.perf_counter()-t_total:.1f}s")

    # Flush remaining buffer
    removed += _flush_buffer()
    bar.close()

    # ── Final report ──────────────────────────────────────────────────
    total_elapsed = time.perf_counter() - t_total
    _log(f"DONE: {removed} removed / {n_edges} total in {total_elapsed:.1f}s")
    _log(f"  precheck: {precheck_necessary} necessary, {precheck_passed} passed")
    _log(f"  apsp: {apsp_stats.get('apsp_calls', 0)} calls, "
         f"total={apsp_stats.get('apsp_total_s', 0):.2f}s")

    remaining = [e for e in sel if exists[e[0], e[1]]]
    return remaining, removed
