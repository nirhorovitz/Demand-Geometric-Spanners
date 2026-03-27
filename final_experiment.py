"""
Yao-Graph Spanner Experiment
=============================
Compare 6 spanner algorithms on random point sets across 4 stretch values.

Algorithms:
  1. greedy(t)              -- standard greedy t-spanner
  2. DGF(complete, t)       -- Descending Greedy Filter on complete graph
  3. yao(t)                 -- t-Yao graph (k equal cones, nearest nbr per cone)
  4. greedy(sqrt(t)) -> DGF(t)
  5. yao(t) -> DGF(t)       -- Yao graph filtered down
  6. yao(t) -> greedy(t, E=yao_t)  -- Yao-seeded greedy

Algorithms 3, 5, 6 share a single precomputed yao(t) per experiment loop.

Usage (smoke test):   change N = 20, then run python experiment_yao.py
Usage (full run):     set N = 5000, then run python experiment_yao.py
"""

from __future__ import annotations

import re
import subprocess
import sys


def _install(pkg: str) -> None:
    print(f"  Installing {pkg}...", flush=True)
    subprocess.check_call([sys.executable, "-m", "pip", "install", pkg],
                          stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def _bootstrap() -> None:
    """Install missing dependencies before any third-party imports."""
    required = [
        ("numpy",      "numpy>=1.24"),
        ("matplotlib", "matplotlib>=3.5"),
        ("tqdm",       "tqdm>=4.0"),
    ]
    for module, spec in required:
        try:
            __import__(module)
        except ImportError:
            _install(spec)

    # CuPy: only attempt if an NVIDIA GPU is present.
    # Detect CUDA version from nvidia-smi, then install the matching wheel.
    try:
        import cupy  # noqa: F401 — already installed, nothing to do
    except ImportError:
        try:
            out = subprocess.run(
                ["nvidia-smi"], capture_output=True, text=True, timeout=10
            ).stdout
            match = re.search(r"CUDA Version:\s*(\d+)\.\d+", out)
            if match:
                major = int(match.group(1))
                if major >= 12:
                    cupy_pkg = "cupy-cuda12x"
                elif major == 11:
                    cupy_pkg = "cupy-cuda11x"
                elif major == 10:
                    cupy_pkg = "cupy-cuda102"
                else:
                    cupy_pkg = None
                if cupy_pkg:
                    _install(cupy_pkg)
        except Exception:
            pass  # No GPU / nvidia-smi not found — CuPy skipped, numpy fallback used


_bootstrap()

import json
import math
import time
from pathlib import Path
from typing import Any

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from tqdm import tqdm

import scipy.sparse as sp_sparse
from scipy.sparse.csgraph import shortest_path as sp_shortest_path, dijkstra as sp_dijkstra

try:
    import cupy as cp
    _CUDA = True
except ImportError:
    _CUDA = False

# ── Project imports (read-only utilities, no pipeline) ────────────────────────
from core.metrics import (
    _euclidean_distances,
    apsp_add_edge,
    apsp_shortest_path_after_removal,
    _compute_degrees,
    _absolute_weight,
    compute_metrics,
)
from core.mst import mst_weight

# ── Plotting constants ────────────────────────────────────────────────────────
NA_STR = "N/A"
TITLE_FONT_SIZE = 14
HEADER_FONT_SIZE = 10
CELL_FONT_SIZE = 9
TABLE_ROW_HEIGHT_IN = 0.34
GRAPH_ROW_HEIGHT_IN = 4.2
COL_WIDTH_IN = 2.8
TABLE_GRID_COLOR = "#B8B8B8"
TABLE_GRID_LINEWIDTH = 0.8
TABLE_WSPACE = 0.06
TABLE_HSPACE = 0.14

# Table columns (in display order)
METRIC_NAMES = [
    "runtime_ms",
    "edge_count",
    "max_stretch",
    "is_minimal",
    "absolute_weight",
    "relative_weight_to_mst",
    "max_degree",
    "average_degree",
]

# ── Plotting helpers (copied verbatim from compressions/plotting.py) ──────────

def _format_metric_value(val: Any, metric_name: str) -> str:
    """Format a metric value for table display."""
    if val is None:
        return NA_STR
    if metric_name == "max_stretch":
        if isinstance(val, str):
            if val in ("Infinity", "-Infinity", "NaN"):
                return val
            return str(val)
        if isinstance(val, (int, float)):
            try:
                if val == float("inf"):
                    return "Infinity"
                if val == float("-inf"):
                    return "-Infinity"
                if math.isnan(val):
                    return "NaN"
            except (TypeError, OverflowError):
                pass
            if isinstance(val, float):
                return f"{val:.8f}"
            return str(val)
        return str(val)
    if metric_name == "status":
        if isinstance(val, bool):
            return "ok" if val else "failed"
        return str(val)
    if metric_name == "is_minimal":
        if isinstance(val, bool):
            return "yes" if val else "no"
        return str(val)
    if metric_name == "runtime_ms":
        if isinstance(val, (int, float)):
            if val == float("inf"):
                return "Infinity"
            if val == float("-inf"):
                return "-Infinity"
            try:
                if math.isnan(val):
                    return "NaN"
            except (TypeError, OverflowError):
                pass
            return f"{float(val) / 1000.0:.3f}s"
        return str(val)
    if isinstance(val, float):
        if val == float("inf"):
            return "Infinity"
        if val == float("-inf"):
            return "-Infinity"
        try:
            if math.isnan(val):
                return "NaN"
        except (TypeError, OverflowError):
            pass
        return f"{val:.3f}"
    return str(val)


def _humanize_metric_name(metric: str) -> str:
    """Convert internal metric key to human-readable label."""
    special_labels = {
        "runtime_ms": "runtime",
        "relative_weight_to_mst": "weight / MST",
        "is_minimal": "minimal",
        "absolute_weight": "abs weight",
        "max_stretch": "max stretch",
        "edge_count": "edges",
        "max_degree": "max deg",
        "average_degree": "avg deg",
    }
    return special_labels.get(metric, metric.replace("_", " "))


def _draw_table_cell(ax: Any, text: str, *, fontsize: int) -> None:
    """Render one table cell with centered text only (grid drawn separately)."""
    ax.axis("off")
    ax.text(0.5, 0.5, text, ha="center", va="center", fontsize=fontsize, wrap=True)


def _draw_table_grid(fig: Any, axes: np.ndarray, table_last_row: int, n_cols: int) -> None:
    """Draw shared separator lines between table rows/columns (no boxed cells)."""
    table_top = axes[0, 0].get_position().y1
    table_bottom = axes[table_last_row, 0].get_position().y0
    table_left = axes[0, 0].get_position().x0
    table_right = axes[0, n_cols - 1].get_position().x1

    for c in range(n_cols - 1):
        x = (axes[0, c].get_position().x1 + axes[0, c + 1].get_position().x0) / 2.0
        fig.add_artist(
            Line2D(
                [x, x],
                [table_bottom, table_top],
                transform=fig.transFigure,
                color=TABLE_GRID_COLOR,
                linewidth=TABLE_GRID_LINEWIDTH,
                zorder=0,
            )
        )

    for r in range(table_last_row):
        y = (axes[r, 0].get_position().y0 + axes[r + 1, 0].get_position().y1) / 2.0
        fig.add_artist(
            Line2D(
                [table_left, table_right],
                [y, y],
                transform=fig.transFigure,
                color=TABLE_GRID_COLOR,
                linewidth=TABLE_GRID_LINEWIDTH,
                zorder=0,
            )
        )


def _normalize_points(points: list[list[float]]) -> np.ndarray:
    """Normalize points to [0,1]x[0,1] for display."""
    arr = np.asarray(points, dtype=np.float64)
    if arr.size == 0:
        return arr.reshape(0, 2)
    if arr.ndim == 1:
        arr = arr.reshape(-1, 2)
    xmin, xmax = arr[:, 0].min(), arr[:, 0].max()
    ymin, ymax = arr[:, 1].min(), arr[:, 1].max()
    xrange = xmax - xmin if xmax > xmin else 1.0
    yrange = ymax - ymin if ymax > ymin else 1.0
    arr_norm = np.zeros_like(arr)
    arr_norm[:, 0] = (arr[:, 0] - xmin) / xrange
    arr_norm[:, 1] = (arr[:, 1] - ymin) / yrange
    return arr_norm


def _canonicalize_edge(i: int, j: int) -> tuple[int, int]:
    """Canonicalize edge to undirected tuple (sorted endpoints)."""
    return (min(i, j), max(i, j))


def _draw_graph_panel(
    ax: Any,
    points: np.ndarray,
    edges: list[list[int]],
    title: str,
    *,
    common_edges: frozenset[tuple[int, int]] | None = None,
    plot_diff: bool = False,
) -> None:
    """Draw one graph panel: points as scatter, edges as segments."""
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=8, pad=2)

    if points.size == 0:
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        return

    pts_norm = _normalize_points(points.tolist())
    ax.scatter(pts_norm[:, 0], pts_norm[:, 1], s=4, c="black", zorder=2)

    if edges:
        n = len(points)
        for edge in edges:
            if len(edge) < 2:
                continue
            i, j = int(edge[0]), int(edge[1])
            if 0 <= i < n and 0 <= j < n:
                if plot_diff and common_edges is not None:
                    canonical = _canonicalize_edge(i, j)
                    color = "r" if canonical not in common_edges else "b"
                else:
                    color = "b"
                ax.plot(
                    [pts_norm[i, 0], pts_norm[j, 0]],
                    [pts_norm[i, 1], pts_norm[j, 1]],
                    f"{color}-",
                    linewidth=0.5,
                    zorder=1,
                    alpha=0.5,
                )

    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xticks([])
    ax.set_yticks([])


# ── JSON serialization helper ─────────────────────────────────────────────────

def _to_jsonable(obj: Any) -> Any:
    """Recursively convert numpy types and infinity to JSON-safe Python types."""
    if isinstance(obj, dict):
        return {k: _to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_to_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _to_jsonable(obj.tolist())
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        f = float(obj)
        if math.isinf(f) or math.isnan(f):
            return None
        return f
    if isinstance(obj, float):
        if math.isinf(obj) or math.isnan(obj):
            return None
        return obj
    if isinstance(obj, np.bool_):
        return bool(obj)
    return obj


# ── Fast vectorized Floyd-Warshall ────────────────────────────────────────────

def _adj_to_scipy(adj_np: np.ndarray) -> np.ndarray:
    """Convert inf-based adjacency to scipy format (0 = no edge). O(n²)."""
    return np.where(np.isinf(adj_np), 0.0, adj_np)


def _fw(adj_np: np.ndarray) -> np.ndarray:
    """
    Floyd-Warshall APSP.  GPU → CuPy vectorized loop.  CPU → scipy C-level FW.
    """
    if _CUDA:
        d = cp.asarray(adj_np, dtype=np.float64)
        n = d.shape[0]
        tmp = cp.empty((n, n), dtype=np.float64)
        for k in range(n):
            cp.add(d[:, k : k + 1], d[k : k + 1, :], out=tmp)
            cp.minimum(d, tmp, out=d)
        return cp.asnumpy(d)
    return sp_shortest_path(_adj_to_scipy(adj_np), method='FW', directed=False)


def _fw_device(adj_np: np.ndarray):
    """Like _fw but keeps result on the xp device — no GPU→CPU transfer."""
    if _CUDA:
        d = cp.asarray(adj_np, dtype=np.float64)
        n = d.shape[0]
        tmp = cp.empty((n, n), dtype=np.float64)
        for k in range(n):
            cp.add(d[:, k:k+1], d[k:k+1, :], out=tmp)
            cp.minimum(d, tmp, out=d)
        return d  # stays on device
    return sp_shortest_path(_adj_to_scipy(adj_np), method='FW', directed=False)


# ── Adjacency matrix builder (shared by DGF and violation check) ──────────────

def _build_adj_np(edge_list: list, n: int, dist_np: np.ndarray) -> np.ndarray:
    """Build a CPU adjacency matrix from a list of (i, j) edge tuples."""
    adj = np.full((n, n), np.inf, dtype=np.float64)
    np.fill_diagonal(adj, 0.0)
    if edge_list:
        E = np.array(edge_list, dtype=np.int64)
        adj[E[:, 0], E[:, 1]] = dist_np[E[:, 0], E[:, 1]]
        adj[E[:, 1], E[:, 0]] = dist_np[E[:, 0], E[:, 1]]
    return adj


def _build_sparse_adj(edge_list: list, n: int, dist_np: np.ndarray) -> sp_sparse.lil_matrix:
    """Build scipy lil_matrix from edge list. O(m) construction."""
    g = sp_sparse.lil_matrix((n, n), dtype=np.float64)
    for u, v in edge_list:
        w = dist_np[u, v]
        g[u, v] = w
        g[v, u] = w
    return g


def _dijkstra_one(adj, src: int, n: int) -> np.ndarray:
    """Dijkstra from src. Accepts sparse (scipy) or dense (numpy) adjacency."""
    if sp_sparse.issparse(adj):
        return sp_dijkstra(adj.tocsr(), indices=[src]).ravel()
    # Dense: convert to sparse (inf → absent, 0 diag → absent)
    finite = np.isfinite(adj) & (adj > 0)
    csr = sp_sparse.csr_matrix(np.where(finite, adj, 0.0))
    return sp_dijkstra(csr, indices=[src]).ravel()


def _dijkstra_with_pred(adj_np: np.ndarray, src: int) -> tuple[np.ndarray, np.ndarray]:
    """Dijkstra from src returning (distances, predecessors). Dense adj input."""
    finite = np.isfinite(adj_np) & (adj_np > 0)
    csr = sp_sparse.csr_matrix(np.where(finite, adj_np, 0.0))
    d, pred = sp_dijkstra(csr, indices=[src], return_predecessors=True)
    return d.ravel(), pred[0]


def _trace_and_update_uses(pred: np.ndarray, r: int, s: int, uses: dict) -> None:
    """Trace SP from r to s via predecessor array. Add (r,s) to uses for each edge on path."""
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


# ── Yao graph ─────────────────────────────────────────────────────────────────

def yao_k_for_t(t: float) -> int:
    """k = ceil(π / arcsin((t−1)/(2t))), minimum 7. No upper cap.

    NOTE: for very small t (e.g. t=1.001), k can be very large (≥ n−1),
    causing the Yao graph to degenerate toward the complete graph.
    """
    val = (t - 1.0) / (2.0 * t)
    # arcsin domain: val must be in (-1, 1). For t > 1 this is always positive < 0.5.
    return max(7, math.ceil(math.pi / math.asin(val)))


def yao_graph(points: np.ndarray, k: int) -> np.ndarray:
    """
    Build undirected Yao_k graph (vectorized).

    For each of k equal-angle cones, find every point's nearest neighbour in
    that cone via vectorized argmin — no Python loop over points.
    Return symmetric closure as sorted (i<j) edge array.

    NOTE: when k >= n-1 (e.g. t=1.001 -> large k, n=5000), the Yao graph
    degenerates to the complete graph; reported as-is.
    """
    n = points.shape[0]
    xp = cp if _CUDA else np
    cone_angle = 2.0 * math.pi / k

    pts = xp.asarray(points)
    x, y = pts[:, 0], pts[:, 1]
    # angle_mat[i,j] = angle of vector from i to j
    dx = x[np.newaxis, :] - x[:, np.newaxis]   # (n, n)
    dy = y[np.newaxis, :] - y[:, np.newaxis]   # (n, n)
    angle_mat = xp.arctan2(dy, dx) % (2.0 * math.pi)

    dist_mat = xp.asarray(_euclidean_distances(points))
    diag = xp.arange(n)
    dist_mat[diag, diag] = xp.inf  # exclude self-edges

    edges: set[tuple[int, int]] = set()
    for c in tqdm(range(k), desc="yao_graph", unit="cone", leave=False):
        lo = c * cone_angle
        hi = (c + 1) * cone_angle
        in_cone = (angle_mat >= lo) & (angle_mat < hi)       # (n, n) bool
        masked = xp.where(in_cone, dist_mat, xp.inf)         # (n, n)
        best_j = xp.argmin(masked, axis=1)                   # (n,)
        best_d = masked[xp.arange(n), best_j]                # (n,)
        valid = ~xp.isinf(best_d)

        valid_np = cp.asnumpy(valid) if _CUDA else np.asarray(valid)
        best_j_np = cp.asnumpy(best_j) if _CUDA else np.asarray(best_j)

        for i in np.where(valid_np)[0]:
            j = int(best_j_np[i])
            edges.add((min(i, j), max(i, j)))

    if not edges:
        return np.empty((0, 2), dtype=np.int64)
    return np.array(sorted(edges), dtype=np.int64)


# ── DGF (clean single-pass) ───────────────────────────────────────────────────

def _check_all_pairs_violation(
    edges_without: list[tuple[int, int]],
    dist_np: np.ndarray,
    w: np.ndarray,   # always ones; kept for API compatibility
    n: int,
    t: float,
) -> bool:
    """
    Return True iff any pair (r, s) violates the t-spanner property in the
    graph defined by edges_without:  sp[r,s] > t * dist[r,s]  OR  sp[r,s] == inf

    Vectorized: builds adjacency with fancy indexing, runs _fw (GPU if available),
    checks the full upper triangle in one operation — no Python pair loop.
    """
    if n < 2:
        return False

    # Build adjacency matrix with fast fancy indexing (always on CPU)
    adj = _build_adj_np(edges_without, n, dist_np)

    # Floyd-Warshall (vectorized, on GPU if available) → numpy result
    sp = _fw(adj)

    # Single vectorized check over upper triangle (w is all-ones)
    r_idx, s_idx = np.triu_indices(n, k=1)
    sp_upper = sp[r_idx, s_idx]
    dist_upper = dist_np[r_idx, s_idx]
    return bool(np.any(np.isinf(sp_upper) | (sp_upper > t * dist_upper)))


def dgf_one_pass(
    edges: list[tuple[int, int]],
    dist_np: np.ndarray,
    n: int,
    t: float,
) -> tuple[list[tuple[int, int]], int]:
    """
    Single-pass Descending Greedy Filter (DGF) — optimized.

    Three speedups over the naive approach:
    1. Batch pre-filter: one FW call removes all edges not on any shortest path.
    2. Incremental adjacency: O(1) update per edge instead of rebuilding.
    3. Dijkstra pre-check: O(n^2) necessity test before O(n^3) FW.

    Returns:
        (remaining_edges, n_removed)
    """
    if len(edges) <= 1:
        return list(edges), 0

    xp = cp if _CUDA else np
    sel = sorted(edges, key=lambda e: -dist_np[e[0], e[1]])

    # Dual-mode: sparse Dijkstra for small edge sets, dense FW for large ones
    use_sparse = len(sel) < 4 * n

    # ── Phase 1: batch pre-filter (one APSP call) ────────────────────────
    if use_sparse:
        lil = _build_sparse_adj(sel, n, dist_np)
        sp_init_np = sp_shortest_path(lil.tocsr(), directed=False)
    else:
        adj = _build_adj_np(sel, n, dist_np)
        sp_init = _fw_device(adj)

    # Vectorized: keep only edges where sp[u,v] >= w (on a shortest path)
    E_arr = np.array(sel, dtype=np.int64)
    if use_sparse:
        sp_vals = sp_init_np[E_arr[:, 0], E_arr[:, 1]]
    else:
        sp_vals = sp_init[E_arr[:, 0], E_arr[:, 1]]
        if _CUDA:
            sp_vals = cp.asnumpy(sp_vals)
    w_vals = dist_np[E_arr[:, 0], E_arr[:, 1]]
    keep_mask = sp_vals >= w_vals  # sp == w means on shortest path
    n_batch_removed = int(np.sum(~keep_mask))

    current_set = set()
    for idx in np.where(keep_mask)[0]:
        e = sel[idx]
        current_set.add((min(e[0], e[1]), max(e[0], e[1])))
    sel_filtered = [e for e in sel if (min(e[0], e[1]), max(e[0], e[1])) in current_set]

    if n_batch_removed > 0:
        tqdm.write(f"    batch pre-filter: removed {n_batch_removed}/{len(sel)} edges, {len(sel_filtered)} remain")

    # ── Phase 2: one-by-one with Dijkstra pre-check ─────────────────────
    # Precompute upper-triangle indices and dist — reused every iteration
    r_np, s_np = np.triu_indices(n, k=1)
    dist_upper_np = dist_np[r_np, s_np]

    def _violation_np(sp_mat):
        """Check violation using numpy array (works for both paths)."""
        sp_up = sp_mat[r_np, s_np]
        return bool(np.any(np.isinf(sp_up) | (sp_up > t * dist_upper_np)))

    if use_sparse:
        # ── Sparse path: scipy Dijkstra (C-level) ───────────────────────
        lil = _build_sparse_adj(sel_filtered, n, dist_np)
        removed = 0

        for edge in tqdm(sel_filtered, desc="dgf(sparse)", unit="edge", leave=False):
            u, v = edge[0], edge[1]
            w = dist_np[u, v]

            # Tentatively remove edge
            lil[u, v] = 0
            lil[v, u] = 0

            # Quick necessity check: Dijkstra from u, O(E log n) in C
            d_from_u = sp_dijkstra(lil.tocsr(), indices=[u]).ravel()
            if d_from_u[v] > t * w:
                lil[u, v] = w
                lil[v, u] = w
                continue

            # Full APSP check via sparse Dijkstra
            sp_new = sp_shortest_path(lil.tocsr(), directed=False)
            if _violation_np(sp_new):
                lil[u, v] = w
                lil[v, u] = w
            else:
                removed += 1

        # Reconstruct remaining edges from sparse matrix
        remaining = []
        for e in sel_filtered:
            u, v = e[0], e[1]
            if lil[u, v] != 0:
                remaining.append(e)
    else:
        # ── Dense path: edge-dependency tracking ─────────────────────────
        adj = _build_adj_np(sel_filtered, n, dist_np)
        removed = 0

        # Initialize uses via APSP with predecessor tracking.
        # For K_n (complete graph), this reduces to uses[(u,v)] = {(u,v)}
        # since direct edges are always shortest (triangle inequality).
        # For non-complete inputs (e.g. greedy output), pairs without
        # direct edges route through multi-hop paths that must be tracked.
        n_edges = len(sel_filtered)
        is_complete = (n_edges == n * (n - 1) // 2)

        uses: dict[tuple[int, int], set[tuple[int, int]]] = {}
        for e in sel_filtered:
            key = (min(e[0], e[1]), max(e[0], e[1]))
            uses[key] = set()

        if is_complete:
            # K_n Euclidean: every pair's SP is its direct edge
            for e in sel_filtered:
                key = (min(e[0], e[1]), max(e[0], e[1]))
                uses[key].add(key)
        else:
            # General graph: run APSP with predecessors, trace all pairs
            finite = np.isfinite(adj) & (adj > 0)
            csr = sp_sparse.csr_matrix(np.where(finite, adj, 0.0))
            _, predecessors = sp_dijkstra(csr, return_predecessors=True)
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

        for edge in tqdm(sel_filtered, desc="dgf", unit="edge", leave=False):
            u, v = edge[0], edge[1]
            w = dist_np[u, v]
            key = (min(u, v), max(u, v))

            pairs = uses.get(key, set())
            if not pairs:
                # Edge not on any current SP — removable instantly
                adj[u, v] = np.inf
                adj[v, u] = np.inf
                removed += 1
                continue

            # Tentatively remove edge
            adj[u, v] = np.inf
            adj[v, u] = np.inf

            # Check each dependent pair via targeted Dijkstra
            necessary = False
            dijkstra_cache: dict[int, tuple[np.ndarray, np.ndarray]] = {}

            for (r, s) in pairs:
                src = r
                if src not in dijkstra_cache:
                    dijkstra_cache[src] = _dijkstra_with_pred(adj, src)
                d_src = dijkstra_cache[src][0]
                if d_src[s] > t * dist_np[r, s]:
                    necessary = True
                    break

            if necessary:
                adj[u, v] = w
                adj[v, u] = w
                continue

            # Removable — update uses for re-routed pairs
            removed += 1
            for (r, s) in pairs:
                if r in dijkstra_cache:
                    pred = dijkstra_cache[r][1]
                    _trace_and_update_uses(pred, r, s, uses)
                elif s in dijkstra_cache:
                    pred = dijkstra_cache[s][1]
                    _trace_and_update_uses(pred, s, r, uses)
                # else: pair was checked via a shared source — its route
                # will be stale in uses but that's harmless (extra checks)
            uses.pop(key, None)

        remaining = []
        for e in sel_filtered:
            u, v = e[0], e[1]
            if not np.isinf(adj[u, v]):
                remaining.append(e)

    total_removed = n_batch_removed + removed
    return remaining, total_removed


# ── Minimality check ──────────────────────────────────────────────────────────

def check_minimality(
    edges: np.ndarray,
    dist: np.ndarray,
    n: int,
    t: float,
) -> bool:
    """
    E is minimal iff dgf_one_pass(E) removes 0 edges.
    (Every edge in E is necessary: removing it breaks some pair's t-path.)
    """
    if len(edges) == 0:
        return True
    edge_list = [(int(e[0]), int(e[1])) for e in edges]
    _, removed = dgf_one_pass(edge_list, dist, n, t)
    return removed == 0


# ── GPU-aware greedy ──────────────────────────────────────────────────────────

def _greedy_local(points: np.ndarray, t: float, E_input=None) -> np.ndarray:
    """
    GPU-aware greedy t-spanner. Replaces greedy_run for final_experiment calls.

    Keeps sp on device; APSP update is one vectorized broadcast per added edge
    instead of the Python double loop in apsp_add_edge.
    """
    from algorithms.base import resolve_candidates
    n = points.shape[0]
    if n <= 1:
        return np.empty((0, 2), dtype=np.int64)

    xp = cp if _CUDA else np
    dist_np = _euclidean_distances(points)
    candidates = resolve_candidates(points, E_input)
    if candidates.shape[0] == 0:
        return np.empty((0, 2), dtype=np.int64)

    lengths = dist_np[candidates[:, 0], candidates[:, 1]]
    candidates = candidates[np.argsort(lengths)]

    selected = []
    sp = None  # APSP on device, initialised on first edge

    for idx in tqdm(range(len(candidates)), desc="greedy", unit="edge", leave=False):
        i, j = int(candidates[idx, 0]), int(candidates[idx, 1])
        w = float(dist_np[i, j])

        if sp is None:
            delta = np.inf
        else:
            val = sp[i, j]
            delta = float(val.get() if _CUDA else val)

        if np.isinf(delta) or delta > t * w:
            selected.append((i, j))
            if sp is None:
                sp = xp.full((n, n), xp.inf, dtype=np.float64)
                sp[xp.arange(n), xp.arange(n)] = 0.0
            # Vectorized APSP update (replaces apsp_add_edge Python double loop)
            cand = xp.minimum(
                sp[:, i:i+1] + w + sp[j:j+1, :],
                sp[:, j:j+1] + w + sp[i:i+1, :]
            )
            xp.minimum(sp, cand, out=sp)

    return np.array(selected, dtype=np.int64) if selected else np.empty((0, 2), dtype=np.int64)


# ── Six algorithm wrappers ────────────────────────────────────────────────────


def algo_greedy(points: np.ndarray, t: float) -> np.ndarray:
    """Algorithm 1: standard greedy t-spanner."""
    return _greedy_local(points, t)


def algo_dgf(points: np.ndarray, t: float, dist: np.ndarray) -> np.ndarray:
    """Algorithm 2: DGF on the complete graph (all n(n-1)/2 pairs)."""
    n = points.shape[0]
    all_pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    result, _ = dgf_one_pass(all_pairs, dist, n, t)
    return np.array(result, dtype=np.int64) if result else np.empty((0, 2), dtype=np.int64)


def algo_yao(yao_t_edges: np.ndarray) -> np.ndarray:
    """Algorithm 3: just the t-Yao graph (precomputed)."""
    return yao_t_edges


def algo_sqrt_greedy_dgf(points: np.ndarray, t: float, dist: np.ndarray) -> np.ndarray:
    """Algorithm 4: greedy(√t) -> DGF(t)."""
    stage1 = _greedy_local(points, math.sqrt(t))
    edge_list = [(int(e[0]), int(e[1])) for e in stage1]
    n = points.shape[0]
    result, _ = dgf_one_pass(edge_list, dist, n, t)
    return np.array(result, dtype=np.int64) if result else np.empty((0, 2), dtype=np.int64)


def algo_yao_dgf(
    yao_t_edges: np.ndarray,
    dist: np.ndarray,
    n: int,
    t: float,
) -> np.ndarray:
    """Algorithm 5: t-Yao -> DGF(t). Reuses precomputed yao_t_edges."""
    edge_list = [(int(e[0]), int(e[1])) for e in yao_t_edges]
    result, _ = dgf_one_pass(edge_list, dist, n, t)
    return np.array(result, dtype=np.int64) if result else np.empty((0, 2), dtype=np.int64)


def algo_yao_greedy_t(
    yao_t_edges: np.ndarray,
    points: np.ndarray,
    t: float,
) -> np.ndarray:
    """Algorithm 6: t-Yao -> greedy(t, E=yao_t). Reuses precomputed yao_t_edges."""
    return _greedy_local(points, t, E_input=yao_t_edges)


# ── Plot helper ───────────────────────────────────────────────────────────────

def _plot_table(
    out_dir: Path,
    points: np.ndarray,
    algo_results: dict[str, dict[str, Any]],
    t: float,
) -> None:
    """
    Render summary_table.png for one experiment (one t-value).

    Layout: rows = algorithms, cols = metrics (transposed from metric-per-row).
    No graph panels.
    """
    algo_names = list(algo_results.keys())
    n_metrics = len(METRIC_NAMES)
    n_algos = len(algo_names)

    # rows: header + one row per algorithm
    # cols: algorithm-name col + one col per metric
    n_rows = n_algos + 1
    n_cols = n_metrics + 1

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(max(6, n_cols * COL_WIDTH_IN), n_rows * TABLE_ROW_HEIGHT_IN),
        squeeze=False,
    )
    fig.suptitle(
        f"Yao Experiment | n={points.shape[0]} | t={t}",
        fontsize=TITLE_FONT_SIZE,
        y=1.02,
    )

    # Header row: blank corner, then one metric per column
    _draw_table_cell(axes[0, 0], "algorithm", fontsize=HEADER_FONT_SIZE)
    for j, metric in enumerate(METRIC_NAMES):
        _draw_table_cell(axes[0, j + 1], _humanize_metric_name(metric), fontsize=HEADER_FONT_SIZE)

    # Algorithm rows
    for i, name in enumerate(algo_names):
        r = i + 1
        _draw_table_cell(axes[r, 0], name, fontsize=HEADER_FONT_SIZE)
        for j, metric in enumerate(METRIC_NAMES):
            val = algo_results[name].get(metric)
            _draw_table_cell(axes[r, j + 1], _format_metric_value(val, metric), fontsize=CELL_FONT_SIZE)

    _draw_table_grid(fig, axes, table_last_row=n_algos, n_cols=n_cols)
    fig.subplots_adjust(top=0.94, hspace=TABLE_HSPACE, wspace=TABLE_WSPACE)

    out_path = out_dir / "summary_table.png"
    fig.savefig(out_path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot saved: {out_path}")


# ── Main experiment loop ──────────────────────────────────────────────────────

T_VALUES = [1.5, 1.1, 1.01, 1.001]
N = 5000
BASE_SEED = 42


def _print_table(algo_results: dict[str, Any], t: float) -> None:
    """Print the results table to the console (rows=algorithms, cols=metrics)."""
    algo_names = list(algo_results.keys())
    headers = ["algorithm"] + [_humanize_metric_name(m) for m in METRIC_NAMES]

    # Build cell strings
    rows = []
    for name in algo_names:
        row = [name]
        for metric in METRIC_NAMES:
            row.append(_format_metric_value(algo_results[name].get(metric), metric))
        rows.append(row)

    # Column widths
    col_widths = [max(len(headers[c]), *(len(rows[r][c]) for r in range(len(rows))))
                  for c in range(len(headers))]

    sep = "+-" + "-+-".join("-" * w for w in col_widths) + "-+"
    def fmt_row(cells):
        return "| " + " | ".join(cells[c].ljust(col_widths[c]) for c in range(len(cells))) + " |"

    print(f"\n  Results for t={t}:")
    print(f"  {sep}")
    print(f"  {fmt_row(headers)}")
    print(f"  {sep}")
    for row in rows:
        print(f"  {fmt_row(row)}")
    print(f"  {sep}\n")


def _save_results(out: Path, algo_results: dict[str, Any]) -> None:
    """Write results.json atomically (temp file + rename) so partial runs are safe."""
    tmp = out / "results.json.tmp"
    tmp.write_text(json.dumps(_to_jsonable(algo_results), indent=2))
    tmp.rename(out / "results.json")


def main() -> None:
    BASE_OUT = Path(f"results_final_experiment_n={N}")
    accel = "GPU (CuPy)" if _CUDA else "CPU (NumPy)"
    print(f"Acceleration: {accel}")
    print(f"Yao experiment: n={N}, t-values={T_VALUES}")
    print(f"Output directory: {BASE_OUT}\n")

    t0_total = time.perf_counter()

    for idx, t in enumerate(T_VALUES):
        print(f"=== t = {t} (seed {BASE_SEED + idx}) ===")
        t0_t = time.perf_counter()

        rng = np.random.default_rng(BASE_SEED + idx)
        P = rng.uniform(0.0, 1.0, (N, 2))

        out = BASE_OUT / f"t={t}"
        out.mkdir(parents=True, exist_ok=True)
        (out / "points.json").write_text(json.dumps(P.tolist()))

        dist = _euclidean_distances(P)

        # Compute yao(t) once — shared by algorithms 3, 5, and 6
        k_t = yao_k_for_t(t)
        print(f"  Yao k={k_t} for t={t}", flush=True)
        t0_yao = time.perf_counter()
        yao_t = yao_graph(P, k_t)
        print(f"  yao_graph: {len(yao_t)} edges, {(time.perf_counter()-t0_yao)*1000:.0f} ms")

        RUNS = [
            ("greedy",          lambda: algo_greedy(P, t)),
            ("dgf",             lambda: algo_dgf(P, t, dist)),
            ("yao",             lambda: algo_yao(yao_t)),
            ("sqrt_greedy_dgf", lambda: algo_sqrt_greedy_dgf(P, t, dist)),
            ("yao_dgf",         lambda: algo_yao_dgf(yao_t, dist, N, t)),
            ("yao_greedy_t",    lambda: algo_yao_greedy_t(yao_t, P, t)),
        ]

        algo_results: dict[str, dict[str, Any]] = {}
        for name, fn in tqdm(RUNS, desc=f"t={t}", unit="algo"):
            print(f"  Running {name}...", end=" ", flush=True)
            t0 = time.perf_counter()
            edges = fn()
            runtime_ms = (time.perf_counter() - t0) * 1000.0

            metrics = compute_metrics(
                edges, P, None, t, runtime_ms,
                skip_full_stretch_above_n=10_000_000,
            )

            print(f"{runtime_ms:.0f} ms, {metrics['edge_count']} edges", end=" ", flush=True)

            is_min = check_minimality(edges, dist, N, t)
            print(f"minimal={is_min}")

            algo_results[name] = {
                **metrics,
                "is_minimal": is_min,
                "edges": edges.tolist(),
            }

            # Save after every algorithm so an interruption loses at most one result
            _save_results(out, algo_results)

        _print_table(algo_results, t)
        _plot_table(out, P, algo_results, t)

        t_elapsed = time.perf_counter() - t0_t
        print(f"  t={t} total: {t_elapsed:.1f}s\n")

    total_elapsed = time.perf_counter() - t0_total
    print(f"All experiments done. Total time: {total_elapsed:.1f}s")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Yao-graph spanner experiment")
    parser.add_argument("--n", type=int, default=N, help="Number of random points (default: 5000)")
    args = parser.parse_args()
    N = args.n
    main()
