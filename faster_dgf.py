"""
Faster DGF (Descending Greedy Filter) v3
=========================================
Clean reimplementation. Two paths:
  - Dense (edges >= 4n): uses-tracking, multi-source Dijkstra on affected pairs
  - Sparse (edges < 4n): full APSP per edge (parallel)

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
    """Zero out (u,v) and (v,u). O(log degree)."""
    for a, b in [(u, v), (v, u)]:
        idx = _csr_find(csr, a, b)
        if idx >= 0:
            csr.data[idx] = 0.0


def _csr_restore(csr: sp_sparse.csr_matrix, u: int, v: int, w: float) -> None:
    """Restore (u,v) and (v,u). O(log degree)."""
    for a, b in [(u, v), (v, u)]:
        idx = _csr_find(csr, a, b)
        if idx >= 0:
            csr.data[idx] = w


def _clean(csr: sp_sparse.csr_matrix) -> sp_sparse.csr_matrix:
    """Copy + eliminate zeros so sp_dijkstra sees removed edges as absent."""
    c = csr.copy()
    c.eliminate_zeros()
    return c


# ── Uses (dependency tracking) ───────────────────────────────────────────────

def _build_uses(csr: sp_sparse.csr_matrix, n: int, edge_keys: set) -> dict:
    """Build uses: edge -> set of pairs whose SP uses that edge."""
    _, pred = sp_dijkstra(csr, return_predecessors=True, directed=False)
    uses: dict[tuple[int, int], set[tuple[int, int]]] = {}
    for key in edge_keys:
        uses[key] = set()
    for r in range(n):
        for s in range(r + 1, n):
            pair = (r, s)
            cur = s
            while cur != r and cur >= 0:
                prev = int(pred[r, cur])
                if prev < 0:
                    break
                ek = (min(prev, cur), max(prev, cur))
                if ek in uses:
                    uses[ek].add(pair)
                cur = prev
    return uses


def _update_uses(
    csr: sp_sparse.csr_matrix,
    pair_list: list[tuple[int, int]],
    uses: dict,
) -> None:
    """Re-trace paths for affected pairs after a removal."""
    sources: list[int] = []
    seen: set[int] = set()
    for r, _ in pair_list:
        if r not in seen:
            seen.add(r)
            sources.append(r)

    _, pred_rows = sp_dijkstra(csr, indices=sources, return_predecessors=True, directed=False)
    pred_rows = np.atleast_2d(pred_rows)
    smap = {src: idx for idx, src in enumerate(sources)}

    for (r, s) in pair_list:
        pair = (min(r, s), max(r, s))
        pred = pred_rows[smap[r]]
        cur = s
        while cur != r and cur >= 0:
            prev = int(pred[cur])
            if prev < 0:
                break
            ek = (min(prev, cur), max(prev, cur))
            if ek in uses:
                uses[ek].add(pair)
            cur = prev


# ── Parallel APSP ────────────────────────────────────────────────────────────

def _get_n_workers() -> int:
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        return os.cpu_count() or 1


def _parallel_apsp(csr: sp_sparse.csr_matrix, n: int, stats: dict | None = None) -> np.ndarray:
    """Multi-threaded APSP. scipy Dijkstra releases GIL."""
    nw = _get_n_workers()

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

    chunk_size = (n + nw - 1) // nw
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


# ── Violation check ──────────────────────────────────────────────────────────

def _any_violation(sp_mat: np.ndarray, dist_np: np.ndarray, n: int, t: float) -> bool:
    """True if any pair violates t-spanner property. GPU accelerated."""
    r_idx, s_idx = np.triu_indices(n, k=1)
    if _CUDA:
        sp_g = cp.asarray(sp_mat)
        d_g = cp.asarray(dist_np)
        return bool(cp.any(cp.isinf(sp_g[r_idx, s_idx]) | (sp_g[r_idx, s_idx] > t * d_g[r_idx, s_idx])))
    sp_up = sp_mat[r_idx, s_idx]
    d_up = dist_np[r_idx, s_idx]
    return bool(np.any(np.isinf(sp_up) | (sp_up > t * d_up)))


# ── Main ─────────────────────────────────────────────────────────────────────

def faster_dgf_one_pass(
    edges: list[tuple[int, int]],
    dist_np: np.ndarray,
    n: int,
    t: float,
) -> tuple[list[tuple[int, int]], int]:
    """
    Descending Greedy Filter.

    Dense: uses-tracking, check only affected pairs.
    Sparse: full parallel APSP per candidate removal.
    """
    if len(edges) <= 1:
        return list(edges), 0

    t_total = time.perf_counter()
    sel = sorted(edges, key=lambda e: -dist_np[e[0], e[1]])
    n_edges = len(sel)
    is_dense = n_edges >= 4 * n

    _log(f"n={n}, edges={n_edges}, t={t}, CUDA={_CUDA}, "
         f"native_v2={_NATIVE_V2}, native_v1={_NATIVE_V1}")
    _log(f"density: {'dense' if is_dense else 'sparse'} "
         f"(threshold=4n={4*n}, edges={n_edges})")

    # ── Native Rust fast path ─────────────────────────────────────────
    if _DGF_NATIVE_ENABLED and is_dense:
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

    # Re-check density
    is_dense = len(sel_filtered) >= 4 * n
    mode = "dense" if is_dense else "sparse"
    nw = _get_n_workers()
    _log(f"Phase 2: mode={mode}, edges={len(sel_filtered)}, "
         f"cpu_workers={nw}, CUDA={'yes' if _CUDA else 'no'}")

    # Build working CSR once
    csr = _build_csr(exists, dist_np, n)

    removed = 0

    # Stats tracking
    stats_phase2: dict = {
        "precheck_necessary": 0,      # edge necessary at Dijkstra p->q
        "precheck_passed": 0,          # passed precheck, need full check
        "full_check_necessary": 0,     # necessary after full check
        "full_check_removable": 0,     # removable after full check
        "free_removals": 0,            # removed via empty uses (dense)
        "skipped_already_removed": 0,  # skipped bc already gone
        "time_precheck_s": 0.0,
        "time_fullcheck_s": 0.0,
        "time_uses_update_s": 0.0,
        "time_clean_s": 0.0,
        "apsp_calls": 0,
        "apsp_total_s": 0.0,
    }
    _last_report = [0]  # mutable for closure
    _report_interval = max(len(sel_filtered) // 20, 50)  # report ~20 times

    def _report_progress(edge_idx: int, force: bool = False) -> None:
        if not force and edge_idx - _last_report[0] < _report_interval:
            return
        _last_report[0] = edge_idx
        elapsed = time.perf_counter() - t_total
        s = stats_phase2
        total_checked = s["precheck_necessary"] + s["precheck_passed"]
        _log(f"progress: edge {edge_idx}/{len(sel_filtered)} | "
             f"removed={removed} | elapsed={elapsed:.1f}s")
        _log(f"  precheck: {s['precheck_necessary']} necessary, "
             f"{s['precheck_passed']} passed -> full check")
        _log(f"  fullcheck: {s['full_check_necessary']} necessary, "
             f"{s['full_check_removable']} removable")
        _log(f"  free={s['free_removals']}, skipped={s['skipped_already_removed']}")
        _log(f"  time: precheck={s['time_precheck_s']:.2f}s, "
             f"fullcheck={s['time_fullcheck_s']:.2f}s, "
             f"uses_update={s['time_uses_update_s']:.2f}s, "
             f"clean={s['time_clean_s']:.2f}s")
        if s.get("apsp_calls", 0) > 0:
            avg = s["apsp_total_s"] / s["apsp_calls"]
            _log(f"  apsp: {s['apsp_calls']} calls, "
                 f"total={s['apsp_total_s']:.2f}s, avg={avg*1000:.1f}ms")
        if is_dense and "uses_size_total" in s:
            _log(f"  uses: {s.get('uses_n_edges',0)} edges tracked, "
                 f"total_pairs={s['uses_size_total']}")

    # ══════════════════════════════════════════════════════════════════
    # DENSE PATH
    # ══════════════════════════════════════════════════════════════════
    if is_dense:
        _log("dense: building initial uses map...")
        t0 = time.perf_counter()
        ccsr = _clean(csr)
        edge_keys = set()
        for e in sel_filtered:
            edge_keys.add((min(e[0], e[1]), max(e[0], e[1])))
        uses = _build_uses(ccsr, n, edge_keys)

        # Batch free removal
        free = [k for k, v in uses.items() if not v]
        for u, v in free:
            exists[u, v] = False
            exists[v, u] = False
            _csr_remove(csr, u, v)
            del uses[(u, v)]
            edge_keys.discard((u, v))
        if free:
            removed += len(free)
            n_prefiltered += len(free)
            stats_phase2["free_removals"] += len(free)

        uses_total = sum(len(v) for v in uses.values())
        _log(f"dense: uses init in {time.perf_counter()-t0:.1f}s | "
             f"batch_free={len(free)} | "
             f"edges_tracked={len(uses)} | total_pairs={uses_total}")

        removed_since_rebuild = 0
        rebuild_interval = max(n // 2, 100)
        _log(f"dense: rebuild_interval={rebuild_interval}")

        for edge_idx, edge in enumerate(tqdm(sel_filtered, desc=f"dgf({mode})", unit="edge", leave=False)):
            p, q = edge[0], edge[1]
            w = dist_np[p, q]
            key = (min(p, q), max(p, q))

            if not exists[p, q]:
                stats_phase2["skipped_already_removed"] += 1
                continue

            pairs = uses.get(key, set())
            if not pairs:
                exists[p, q] = False
                exists[q, p] = False
                _csr_remove(csr, p, q)
                uses.pop(key, None)
                edge_keys.discard(key)
                removed += 1
                removed_since_rebuild += 1
                stats_phase2["free_removals"] += 1
                continue

            # Tentatively remove
            exists[p, q] = False
            exists[q, p] = False
            _csr_remove(csr, p, q)

            t0c = time.perf_counter()
            ccsr = _clean(csr)
            stats_phase2["time_clean_s"] += time.perf_counter() - t0c

            # Quick check: Dijkstra p->q
            t0p = time.perf_counter()
            d_pq = sp_dijkstra(ccsr, indices=[p], limit=t * w, directed=False).ravel()
            stats_phase2["time_precheck_s"] += time.perf_counter() - t0p

            if d_pq[q] > t * w:
                exists[p, q] = True
                exists[q, p] = True
                _csr_restore(csr, p, q, w)
                stats_phase2["precheck_necessary"] += 1
                _report_progress(edge_idx)
                continue

            stats_phase2["precheck_passed"] += 1

            # Check affected pairs
            t0f = time.perf_counter()
            pair_list = list(pairs)
            if key in pairs:
                pair_list = [key] + [x for x in pair_list if x != key]

            sources: list[int] = []
            seen: set[int] = set()
            for r, _ in pair_list:
                if r not in seen:
                    seen.add(r)
                    sources.append(r)

            d_rows = sp_dijkstra(ccsr, indices=sources, directed=False)
            d_rows = np.atleast_2d(d_rows)
            smap = {src: idx for idx, src in enumerate(sources)}

            necessary = False
            for (r, s) in pair_list:
                if d_rows[smap[r], s] > t * dist_np[r, s]:
                    necessary = True
                    break
            stats_phase2["time_fullcheck_s"] += time.perf_counter() - t0f

            if necessary:
                exists[p, q] = True
                exists[q, p] = True
                _csr_restore(csr, p, q, w)
                stats_phase2["full_check_necessary"] += 1
                _report_progress(edge_idx)
                continue

            # Removable
            removed += 1
            removed_since_rebuild += 1
            uses.pop(key, None)
            edge_keys.discard(key)
            stats_phase2["full_check_removable"] += 1

            t0u = time.perf_counter()
            _update_uses(ccsr, pair_list, uses)
            stats_phase2["time_uses_update_s"] += time.perf_counter() - t0u

            _report_progress(edge_idx)

            # Periodic rebuild
            if removed_since_rebuild >= rebuild_interval:
                _log(f"dense: rebuilding uses (removed_since={removed_since_rebuild})...")
                t0r = time.perf_counter()
                ccsr = _clean(csr)
                uses = _build_uses(ccsr, n, edge_keys)
                removed_since_rebuild = 0
                nf = 0
                for k in [k for k, v in uses.items() if not v]:
                    u, v = k
                    exists[u, v] = False
                    exists[v, u] = False
                    _csr_remove(csr, u, v)
                    del uses[k]
                    edge_keys.discard(k)
                    nf += 1
                if nf:
                    removed += nf
                    stats_phase2["free_removals"] += nf
                uses_total = sum(len(v) for v in uses.values())
                stats_phase2["uses_size_total"] = uses_total
                stats_phase2["uses_n_edges"] = len(uses)
                _log(f"dense: rebuild done in {time.perf_counter()-t0r:.1f}s | "
                     f"batch_free={nf} | edges={len(uses)} | pairs={uses_total}")

    # ══════════════════════════════════════════════════════════════════
    # SPARSE PATH
    # ══════════════════════════════════════════════════════════════════
    else:
        _log(f"sparse: {nw} CPU workers for parallel APSP, "
             f"n={n} sources per APSP")

        # First APSP to get baseline timing
        t0_first = time.perf_counter()
        first_stats: dict = {}
        _first_test = _parallel_apsp(_clean(csr), n, first_stats)
        del _first_test
        _log(f"sparse: test APSP took {time.perf_counter()-t0_first:.3f}s "
             f"(method={first_stats.get('apsp_method')}, "
             f"workers={first_stats.get('apsp_workers')}, "
             f"chunks={first_stats.get('apsp_chunks')}, "
             f"chunk_sizes={first_stats.get('apsp_chunk_sizes')})")
        if "apsp_chunk_min_s" in first_stats:
            _log(f"sparse: chunk times: "
                 f"min={first_stats['apsp_chunk_min_s']*1000:.1f}ms, "
                 f"max={first_stats['apsp_chunk_max_s']*1000:.1f}ms, "
                 f"avg={first_stats['apsp_chunk_avg_s']*1000:.1f}ms")

        for edge_idx, edge in enumerate(tqdm(sel_filtered, desc=f"dgf({mode})", unit="edge", leave=False)):
            p, q = edge[0], edge[1]
            w = dist_np[p, q]

            if not exists[p, q]:
                stats_phase2["skipped_already_removed"] += 1
                continue

            # Tentatively remove
            exists[p, q] = False
            exists[q, p] = False
            _csr_remove(csr, p, q)

            t0c = time.perf_counter()
            ccsr = _clean(csr)
            stats_phase2["time_clean_s"] += time.perf_counter() - t0c

            # Quick check: Dijkstra p->q
            t0p = time.perf_counter()
            d_pq = sp_dijkstra(ccsr, indices=[p], limit=t * w, directed=False).ravel()
            stats_phase2["time_precheck_s"] += time.perf_counter() - t0p

            if d_pq[q] > t * w:
                exists[p, q] = True
                exists[q, p] = True
                _csr_restore(csr, p, q, w)
                stats_phase2["precheck_necessary"] += 1
                _report_progress(edge_idx)
                continue

            stats_phase2["precheck_passed"] += 1

            # Full APSP check (parallel)
            t0f = time.perf_counter()
            apsp_stats: dict = {}
            sp_all = _parallel_apsp(ccsr, n, apsp_stats)
            stats_phase2["apsp_calls"] = stats_phase2.get("apsp_calls", 0) + 1
            stats_phase2["apsp_total_s"] = stats_phase2.get("apsp_total_s", 0.0) + apsp_stats.get("apsp_last_s", 0.0)

            if _any_violation(sp_all, dist_np, n, t):
                exists[p, q] = True
                exists[q, p] = True
                _csr_restore(csr, p, q, w)
                stats_phase2["full_check_necessary"] += 1
            else:
                removed += 1
                stats_phase2["full_check_removable"] += 1
            stats_phase2["time_fullcheck_s"] += time.perf_counter() - t0f

            _report_progress(edge_idx)

    # ── Final report ──────────────────────────────────────────────────
    total_elapsed = time.perf_counter() - t_total
    _log(f"DONE: {n_prefiltered + removed} total removed "
         f"({n_prefiltered} prefilter + {removed} phase2) in {total_elapsed:.1f}s")
    _report_progress(len(sel_filtered), force=True)

    remaining = [e for e in sel if exists[e[0], e[1]]]
    return remaining, n_prefiltered + removed
