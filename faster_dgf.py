"""
Faster DGF (Descending Greedy Filter) v5
=========================================
Per-edge GPU APSP approach: for each candidate edge, do a quick Dijkstra
precheck, then full GPU Floyd-Warshall APSP + violation check.

Interface:
    faster_dgf_one_pass(edges, dist_np, n, t) -> (remaining_edges, n_removed)
"""

from __future__ import annotations

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



# ── Main ─────────────────────────────────────────────────────────────────────

def faster_dgf_one_pass(
    edges: list[tuple[int, int]],
    dist_np: np.ndarray,
    n: int,
    t: float,
) -> tuple[list[tuple[int, int]], int]:
    """
    Descending Greedy Filter.

    Per-edge: Dijkstra precheck, then full APSP + violation check.
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

    # ── exists matrix ─────────────────────────────────────────────────
    exists = np.zeros((n, n), dtype=bool)
    for u, v in sel:
        exists[u, v] = True
        exists[v, u] = True

    # ── Phase 1: batch pre-filter ─────────────────────────────────────
    is_complete = (n_edges == n * (n - 1) // 2)
    n_prefiltered = 0

    if is_complete:
        _log(f"K_n detected ({n_edges} edges): skipping Phase 1")
        sel_filtered = list(sel)
    else:
        _log("Phase 1: computing APSP for batch pre-filter...")
        t0 = time.perf_counter()
        csr = _build_csr(exists, dist_np, n)
        stats: dict = {}
        sp_init = _parallel_apsp(csr, n, stats)
        _log(f"Phase 1 APSP: {time.perf_counter()-t0:.3f}s "
             f"(method={stats.get('apsp_method')}, workers={stats.get('apsp_workers')})")

        E = np.array(sel, dtype=np.int64)
        sp_vals = sp_init[E[:, 0], E[:, 1]]
        w_vals = dist_np[E[:, 0], E[:, 1]]
        keep = sp_vals >= w_vals - 1e-12

        n_prefiltered = int(np.sum(~keep))
        if n_prefiltered > 0:
            for idx in np.where(~keep)[0]:
                u, v = sel[idx]
                exists[u, v] = False
                exists[v, u] = False
        _log(f"Phase 1: removed {n_prefiltered}/{n_edges}, "
             f"{int(np.sum(keep))} remain")

        sel_filtered = [sel[i] for i in np.where(keep)[0]]

    # ── Phase 2: precheck + GPU APSP per edge ─────────────────────────
    _log(f"Phase 2: precheck + APSP per edge, edges={len(sel_filtered)}, "
         f"CUDA={'yes' if _CUDA else 'no'}")

    # Build working CSR once
    csr = _build_csr(exists, dist_np, n)

    removed = 0

    # Stats tracking
    stats_phase2: dict = {
        "precheck_necessary": 0,      # edge necessary at Dijkstra p->q
        "precheck_passed": 0,          # passed precheck, need full APSP
        "full_check_necessary": 0,     # necessary after full APSP check
        "full_check_removable": 0,     # removable after full APSP check
    }
    apsp_stats: dict = {}
    _last_report = [0]
    _report_interval = max(len(sel_filtered) // 20, 50)

    def _report_progress(edge_idx: int, force: bool = False) -> None:
        if not force and edge_idx - _last_report[0] < _report_interval:
            return
        _last_report[0] = edge_idx
        elapsed = time.perf_counter() - t_total
        s = stats_phase2
        _log(f"progress: edge {edge_idx}/{len(sel_filtered)} | "
             f"removed={removed} | elapsed={elapsed:.1f}s")
        _log(f"  precheck: {s['precheck_necessary']} necessary, "
             f"{s['precheck_passed']} passed -> full APSP")
        _log(f"  fullcheck: {s['full_check_necessary']} necessary, "
             f"{s['full_check_removable']} removable")
        if apsp_stats.get("apsp_calls", 0) > 0:
            avg = apsp_stats["apsp_total_s"] / apsp_stats["apsp_calls"]
            _log(f"  apsp: {apsp_stats['apsp_calls']} calls, "
                 f"total={apsp_stats['apsp_total_s']:.2f}s, avg={avg*1000:.1f}ms")

    # ══════════════════════════════════════════════════════════════════
    # Per-edge: precheck + full APSP + violation check
    # ══════════════════════════════════════════════════════════════════
    for edge_idx, edge in enumerate(tqdm(sel_filtered, desc="dgf(apsp)", unit="edge", leave=False)):
        p, q = edge[0], edge[1]
        w = dist_np[p, q]

        if not exists[p, q]:
            continue

        # 1. Tentatively remove
        _csr_remove(csr, p, q)

        # 2. Precheck: Dijkstra from p with limit = t*w
        d_pq = sp_dijkstra(csr, indices=[p], limit=t * w, directed=False).ravel()

        if d_pq[q] > t * w:
            # Edge is necessary — restore immediately
            _csr_restore(csr, p, q, w)
            stats_phase2["precheck_necessary"] += 1
            _report_progress(edge_idx)
            continue

        stats_phase2["precheck_passed"] += 1

        # 3. Full APSP on current CSR
        sp_matrix = _parallel_apsp(csr, n, apsp_stats)

        # 4. Violation check: any sp[i,j] > t * dist[i,j]?
        if _any_violation(sp_matrix, dist_np, t, n):
            # Violation found — edge is necessary, restore
            _csr_restore(csr, p, q, w)
            stats_phase2["full_check_necessary"] += 1
        else:
            # No violation — edge permanently removed
            exists[p, q] = False
            exists[q, p] = False
            removed += 1
            stats_phase2["full_check_removable"] += 1

        _report_progress(edge_idx)

    # ── Final report ──────────────────────────────────────────────────
    total_elapsed = time.perf_counter() - t_total
    total_removed = n_prefiltered + removed
    _log(f"DONE: {total_removed} total removed "
         f"({n_prefiltered} prefilter + {removed} phase2) in {total_elapsed:.1f}s")
    _report_progress(len(sel_filtered), force=True)

    remaining = [e for e in sel if exists[e[0], e[1]]]
    return remaining, total_removed
