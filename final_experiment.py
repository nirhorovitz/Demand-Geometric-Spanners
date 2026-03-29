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
import shutil
import subprocess
import sys
import os


def _install(pkg: str) -> None:
    print(f"  Installing {pkg}...", flush=True)
    subprocess.check_call([sys.executable, "-m", "pip", "install", pkg],
                          stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def _ensure_rust() -> bool:
    """Install Rust via rustup if cargo is not on PATH. Returns True if available."""
    if shutil.which("cargo"):
        return True
    # Try common rustup install location
    cargo_home = os.path.expanduser("~/.cargo/bin/cargo")
    if os.path.isfile(cargo_home):
        os.environ["PATH"] = os.path.expanduser("~/.cargo/bin") + os.pathsep + os.environ.get("PATH", "")
        return True
    print("  Installing Rust via rustup (needed once for native DGF)...", flush=True)
    try:
        subprocess.check_call(
            ["sh", "-c", "curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=300,
        )
        os.environ["PATH"] = os.path.expanduser("~/.cargo/bin") + os.pathsep + os.environ.get("PATH", "")
        return shutil.which("cargo") is not None
    except Exception as exc:
        print(f"  Rust install failed: {exc}", flush=True)
        return False


def _build_native_dgf() -> bool:
    """Build native_dgf extension in-place using maturin. Returns True on success."""
    native_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "native_dgf")
    if not os.path.isdir(native_dir):
        return False
    try:
        __import__("native_dgf")
        return True  # already importable
    except ImportError:
        pass
    if not _ensure_rust():
        print("  Skipping native DGF build (no Rust toolchain)", flush=True)
        return False
    # Ensure maturin is installed
    try:
        __import__("maturin")
    except ImportError:
        print("  Installing maturin...", flush=True)
        try:
            subprocess.check_call(
                [sys.executable, "-m", "pip", "install", "maturin>=1.7,<2.0"],
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=120,
            )
        except Exception as exc:
            print(f"  maturin install failed: {exc}", flush=True)
            return False
    print("  Building native_dgf (Rust extension)...", flush=True)
    try:
        subprocess.check_call(
            [sys.executable, "-m", "maturin", "develop", "--release", "--manifest-path",
             os.path.join(native_dir, "Cargo.toml")],
            timeout=600,
        )
        return True
    except Exception as exc:
        print(f"  native_dgf build failed: {exc}", flush=True)
        return False


def _bootstrap() -> None:
    """Install missing dependencies before any third-party imports."""
    required = [
        ("numpy",      "numpy>=1.24"),
        ("scipy",      "scipy>=1.9"),
        ("matplotlib", "matplotlib>=3.5"),
        ("tqdm",       "tqdm>=4.0"),
    ]
    for module, spec in required:
        try:
            __import__(module)
        except ImportError:
            _install(spec)

    # CuPy: only attempt if an NVIDIA GPU is present.
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
            pass  # No GPU / nvidia-smi not found — CuPy skipped

    # Native Rust DGF extension — build if not already available
    _build_native_dgf()


_bootstrap()

import json
import math
import time
from pathlib import Path
from typing import Any, Optional

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from tqdm import tqdm

import scipy.sparse as sp_sparse
from scipy.sparse.csgraph import shortest_path as sp_shortest_path, dijkstra as sp_dijkstra

try:
    import cupy as cp
    cp.array([1.0])  # verify GPU is actually usable
    _CUDA = True
except Exception:
    _CUDA = False

try:
    from native_dgf import (
        trace_path_edge_ids as native_trace_path_edge_ids,
        sort_pairs_with_key_first as native_sort_pairs_with_key_first,
        dgf_one_pass_native as native_dgf_one_pass,
    )
    _NATIVE_DGF_AVAILABLE = True
except Exception:
    _NATIVE_DGF_AVAILABLE = False

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

    def report(self, prefix: str = "dgf-profile") -> None:
        if not self.enabled:
            return
        total = sum(self._marks.values())
        print(f"  [{prefix}] total={total:.3f}s")
        for k in sorted(self._marks.keys()):
            v = self._marks[k]
            pct = (100.0 * v / total) if total > 0 else 0.0
            print(f"  [{prefix}] {k}={v:.3f}s ({pct:.1f}%)")

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


def _dijkstra_with_pred_sparse(sparse_adj: sp_sparse.spmatrix, src: int) -> tuple[np.ndarray, np.ndarray]:
    """Dijkstra from src on a sparse adjacency matrix (avoids dense→sparse conversion)."""
    d, pred = sp_dijkstra(sparse_adj.tocsr(), indices=[src], return_predecessors=True)
    return d.ravel(), pred[0]


def _rebuild_uses(sparse_adj: sp_sparse.spmatrix, n: int, existing_keys: set) -> dict:
    """Rebuild uses mapping from scratch via full APSP with predecessors."""
    csr = sparse_adj.tocsr()
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


def _sort_pairs_with_key_first(
    pairs: set[tuple[int, int]],
    key: tuple[int, int],
    *,
    use_native: bool,
) -> list[tuple[int, int]]:
    """Return pair list with `key` first when present, then remaining pairs."""
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


def _trace_and_update_uses_maybe_native(
    pred: np.ndarray,
    r: int,
    s: int,
    uses: dict,
    n: int,
    *,
    use_native: bool,
) -> None:
    """
    Trace SP path and update uses.
    Falls back to Python tracing if native helper is unavailable.
    """
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

    sel = sorted(edges, key=lambda e: -dist_np[e[0], e[1]])

    if _NATIVE_DGF_AVAILABLE and _DGF_NATIVE_ENABLED:
        try:
            return native_dgf_one_pass(sel, dist_np, n, t)
        except Exception as e:
            print(f"    [Native DGF failed: {e}, falling back to Python]")

    xp = cp if _CUDA else np
    profiler = _DgfProfiler(_DGF_PROFILE_ENABLED)
    t_phase = time.perf_counter()
    native_mode = False

    # Dual-mode: sparse Dijkstra for small edge sets, dense FW for large ones
    use_sparse = len(sel) < 4 * n

    # ── Phase 1: batch pre-filter (one APSP call) ────────────────────────
    n_edges = len(sel)
    is_complete = (n_edges == n * (n - 1) // 2)

    if is_complete:
        # K_n APSP skip: for complete graphs, direct edges are always shortest
        # paths (triangle inequality), so Phase 1 cannot remove any edges.
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
        # Build sel_filtered from keep_mask (K_n already set sel_filtered above)
        current_set = set()
        for idx in np.where(keep_mask)[0]:
            e = sel[idx]
            current_set.add((min(e[0], e[1]), max(e[0], e[1])))
        sel_filtered = [e for e in sel if (min(e[0], e[1]), max(e[0], e[1])) in current_set]

        if n_batch_removed > 0:
            tqdm.write(f"    batch pre-filter: removed {n_batch_removed}/{len(sel)} edges, {len(sel_filtered)} remain")
        profiler.add("phase1_prefilter_filtering", time.perf_counter() - t_phase)
        t_phase = time.perf_counter()

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
        # ── Dense path: GPU-accelerated with edge-dependency tracking ─────
        # GPU FW only pays off when n is large enough to amortise kernel launch
        # overhead.  On RTX 3090 + CUDA 12.4, breakeven is around n ≈ 100.
        use_gpu_phase2 = _CUDA and n >= 100

        adj = _build_adj_np(sel_filtered, n, dist_np)
        sparse_adj = _build_sparse_adj(sel_filtered, n, dist_np)
        removed = 0
        rebuild_interval = max(n, 500)

        use_native_helpers = bool(is_complete and _NATIVE_DGF_AVAILABLE and _DGF_NATIVE_ENABLED)
        native_mode = use_native_helpers

        # Dependency map: edge -> set of pairs currently routed through it.
        uses: dict[tuple[int, int], set[tuple[int, int]]] = {}

        if not is_complete:
            for e in sel_filtered:
                key = (min(e[0], e[1]), max(e[0], e[1]))
                uses[key] = set()
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

        shorter_half_start = len(sel_filtered) // 2

        # GPU Phase 2: keep APSP on GPU for fast violation checks
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

        for edge_idx, edge in enumerate(tqdm(sel_filtered, desc="dgf", unit="edge", leave=False)):
            u, v = edge[0], edge[1]
            w = dist_np[u, v]
            key = (min(u, v), max(u, v))

            pairs_explicit = uses.get(key, set())
            if is_complete:
                if pairs_explicit:
                    pair_set = set(pairs_explicit)
                    pair_set.add(key)
                    pair_list = _sort_pairs_with_key_first(
                        pair_set, key, use_native=use_native_helpers
                    )
                else:
                    pair_list = [key]
            else:
                if not pairs_explicit:
                    adj[u, v] = np.inf
                    adj[v, u] = np.inf
                    sparse_adj[u, v] = 0
                    sparse_adj[v, u] = 0
                    if use_gpu_phase2:
                        adj_gpu[u, v] = cp.inf
                        adj_gpu[v, u] = cp.inf
                    removed += 1
                    continue
                pair_list = _sort_pairs_with_key_first(
                    pairs_explicit, key, use_native=use_native_helpers
                )

            # Tentatively remove edge
            adj[u, v] = np.inf
            adj[v, u] = np.inf
            sparse_adj[u, v] = 0
            sparse_adj[v, u] = 0
            if use_gpu_phase2:
                adj_gpu[u, v] = cp.inf
                adj_gpu[v, u] = cp.inf

            # ── GPU-accelerated necessity check ──────────────────────────
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

                # Check which pairs are affected by this edge removal
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
                    uses.pop(key, None)
                    profiler.add("gpu_phase2_free_removal", time.perf_counter() - t0)
                    continue

                # Decide: GPU FW (many sources) vs CPU Dijkstra (few sources)
                # On RTX 3090, GPU FW (O(n³) vectorised) beats CPU Dijkstra
                # from k sources once k is even modestly large.  The GPU
                # also gives us a full APSP update we can cache, so prefer
                # it for anything beyond a handful of sources.
                n_sources = len(set(r for r, _ in pair_list))
                use_gpu_fw = (n_sources > n // 8)

                if use_gpu_fw:
                    # Full GPU FW recompute
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
                        sparse_adj[u, v] = w
                        sparse_adj[v, u] = w
                        adj_gpu[u, v] = w
                        adj_gpu[v, u] = w
                        continue
                    else:
                        sp_gpu = d_gpu
                        removed += 1

                        # Uses update (CPU — need predecessors)
                        t0 = time.perf_counter()
                        csr_current = sparse_adj.tocsr(copy=False)
                        sources_gpu: list[int] = []
                        seen_gpu: set[int] = set()
                        for r, _ in pair_list:
                            if r not in seen_gpu:
                                seen_gpu.add(r)
                                sources_gpu.append(r)
                        _, pred_rows = sp_dijkstra(
                            csr_current, indices=sources_gpu, return_predecessors=True
                        )
                        pred_rows = np.atleast_2d(pred_rows)
                        source_row = {src: idx for idx, src in enumerate(sources_gpu)}
                        for (r, s) in pair_list:
                            pred_src = pred_rows[source_row[r]]
                            _trace_and_update_uses_maybe_native(
                                pred_src, r, s, uses, n, use_native=use_native_helpers
                            )
                        uses.pop(key, None)
                        profiler.add("gpu_phase2_uses_update", time.perf_counter() - t0)

                        if removed > 0 and removed % rebuild_interval == 0:
                            t0 = time.perf_counter()
                            existing_keys = set(uses.keys())
                            uses = _rebuild_uses(sparse_adj, n, existing_keys)
                            profiler.add("dense_uses_rebuild", time.perf_counter() - t0)
                        continue
                else:
                    # Few sources — fall through to CPU Dijkstra path
                    pass

            # ── CPU fallback (non-GPU or few-sources path) ───────────────
            t0 = time.perf_counter()
            csr_current = sparse_adj.tocsr(copy=False)
            profiler.add("dense_csr_build", time.perf_counter() - t0)

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

            t0 = time.perf_counter()
            if split_pred_mode:
                d_rows = sp_dijkstra(csr_current, indices=sources)
            else:
                d_rows, pred_rows = sp_dijkstra(
                    csr_current, indices=sources, return_predecessors=True
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
                sparse_adj[u, v] = w
                sparse_adj[v, u] = w
                if use_gpu_phase2:
                    adj_gpu[u, v] = w
                    adj_gpu[v, u] = w
                continue

            # Removable — update uses for re-routed pairs
            removed += 1
            if use_gpu_phase2:
                sp_gpu = None  # Invalidate GPU APSP cache
            t0 = time.perf_counter()
            if pred_rows is None:
                _, pred_rows = sp_dijkstra(
                    csr_current, indices=sources, return_predecessors=True
                )
            pred_rows = np.atleast_2d(pred_rows)

            for (r, s) in pair_list:
                pred_src = pred_rows[source_row[r]]
                _trace_and_update_uses_maybe_native(
                    pred_src, r, s, uses, n, use_native=use_native_helpers
                )
            uses.pop(key, None)
            profiler.add("dense_uses_update", time.perf_counter() - t0)

            if removed > 0 and removed % rebuild_interval == 0:
                t0 = time.perf_counter()
                existing_keys = set(uses.keys())
                uses = _rebuild_uses(sparse_adj, n, existing_keys)
                profiler.add("dense_uses_rebuild", time.perf_counter() - t0)

        remaining = []
        for e in sel_filtered:
            u, v = e[0], e[1]
            if not np.isinf(adj[u, v]):
                remaining.append(e)

    total_removed = n_batch_removed + removed
    profiler.add("phase2_total", time.perf_counter() - t_phase)
    _gpu_used = not use_sparse and _CUDA and n >= 100
    if _DGF_PROFILE_ENABLED:
        mode = "sparse" if use_sparse else "dense"
        if native_mode:
            mode = f"{mode}-native"
        if _gpu_used:
            mode = f"{mode}-gpu"
        profiler.report(prefix=f"dgf-{mode}")
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

    When E_input is provided (greedy on a t-spanner), uses the modified check:
        add edge iff  δ_{E_curr}(p,q) > √t · λ_{E_input}(p,q)
    where λ_{E_input} is the shortest-path distance in the input graph.
    Otherwise (E_input=None, full graph), uses the standard check:
        add edge iff  δ_{E_curr}(p,q) > t · |pq|
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

    # When E_input is provided, precompute APSP of the input graph
    # and use the modified threshold: √t · λ_{E_input}(p,q)
    sp_input = None
    sqrt_t = math.sqrt(t)
    if E_input is not None:
        adj_input = _build_adj_np(
            [(int(e[0]), int(e[1])) for e in E_input], n, dist_np
        )
        sp_input = _fw(adj_input)

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

        if sp_input is not None:
            threshold = sqrt_t * sp_input[i, j]
        else:
            threshold = t * w

        if np.isinf(delta) or delta > threshold:
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

        out = BASE_OUT / f"t={t}"
        out.mkdir(parents=True, exist_ok=True)

        # --- Resume detection ---
        results_file = out / "results.json"
        points_file = out / "points.json"
        algo_results: dict[str, dict[str, Any]] = {}

        if results_file.exists():
            algo_results = json.loads(results_file.read_text())

        if points_file.exists():
            P = np.array(json.loads(points_file.read_text()))
        else:
            rng = np.random.default_rng(BASE_SEED + idx)
            P = rng.uniform(0.0, 1.0, (N, 2))
            points_file.write_text(json.dumps(P.tolist()))

        dist = _euclidean_distances(P)

        # Compute yao(t) once — shared by algorithms 3, 5, and 6
        yao_t = None  # default; computed below only if needed

        RUNS = [
            ("greedy",          lambda: algo_greedy(P, t)),
            ("dgf",             lambda: algo_dgf(P, t, dist)),
            ("yao",             lambda: algo_yao(yao_t)),
            ("sqrt_greedy_dgf", lambda: algo_sqrt_greedy_dgf(P, t, dist)),
            ("yao_dgf",         lambda: algo_yao_dgf(yao_t, dist, N, t)),
            ("yao_greedy_t",    lambda: algo_yao_greedy_t(yao_t, P, t)),
        ]

        done = set(algo_results.keys())

        # Skip K_n DGF for n >= 1000 — O(n^3.3) is prohibitive
        # if N >= 1000:
        #     if "dgf" not in done:
        #         print(f"  Skipping K_n DGF (n={N} >= 1000, too slow)")
        #     RUNS = [(name, fn) for name, fn in RUNS if name != "dgf"]

        missing_runs = [(name, fn) for name, fn in RUNS if name not in done]

        if not missing_runs:
            print(f"  All algorithms already done, skipping")
            _print_table(algo_results, t)
            t_elapsed = time.perf_counter() - t0_t
            print(f"  t={t} total: {t_elapsed:.1f}s\n")
            continue

        if done:
            print(f"  Resuming: {len(done)}/6 done, running {len(missing_runs)} remaining")

        # Compute yao graph only if any yao-dependent algorithm is missing
        yao_dependent = {"yao", "yao_dgf", "yao_greedy_t"}
        if any(name in yao_dependent for name, _ in missing_runs):
            k_t = yao_k_for_t(t)
            print(f"  Yao k={k_t} for t={t}", flush=True)
            t0_yao = time.perf_counter()
            yao_t = yao_graph(P, k_t)
            print(f"  yao_graph: {len(yao_t)} edges, {(time.perf_counter()-t0_yao)*1000:.0f} ms")

        for name, fn in tqdm(missing_runs, desc=f"t={t}", unit="algo"):
            print(f"  Running {name}...", end=" ", flush=True)
            t0 = time.perf_counter()
            edges = fn()
            runtime_ms = (time.perf_counter() - t0) * 1000.0

            metrics = compute_metrics(
                edges, P, None, t, runtime_ms,
                skip_full_stretch_above_n=10_000_000,
            )

            print(f"{runtime_ms:.0f} ms, {metrics['edge_count']} edges", end=" ", flush=True)

            # DGF-based algorithms already produce minimal spanners by
            # construction — skip the expensive second DGF pass for them.
            dgf_minimal_algos = {"sqrt_greedy_dgf", "yao_dgf", "dgf"}
            if name in dgf_minimal_algos:
                is_min = True
                print(f"minimal={is_min} (by construction)")
            else:
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
