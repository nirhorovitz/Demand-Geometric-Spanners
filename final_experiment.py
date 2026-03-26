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

import json
import math
import time
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from tqdm import tqdm

# ── Project imports (read-only utilities, no pipeline) ────────────────────────
from core.metrics import (
    _euclidean_distances,
    _build_adjacency,
    _floyd_warshall,
    apsp_add_edge,
    apsp_shortest_path_after_removal,
    _compute_degrees,
    _absolute_weight,
    compute_metrics,
)
from core.mst import mst_weight
from algorithms.greedy import run as greedy_run

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
    Build undirected Yao_k graph.

    For each point p and each of k equal-angle cones, add directed edge
    p -> nearest neighbor in that cone (if any). Return symmetric closure
    as sorted (i<j) edge array.

    NOTE: when k >= n-1 (e.g. t=1.001 -> large k, n=5000), the Yao graph
    degenerates to the complete graph; reported as-is.
    """
    n = points.shape[0]
    cone_angle = 2.0 * math.pi / k
    dist = _euclidean_distances(points)
    edges: set[tuple[int, int]] = set()
    for p in tqdm(range(n), desc="yao_graph", unit="pt", leave=False):
        for c in range(k):
            lo = c * cone_angle
            hi = (c + 1) * cone_angle
            best_q, best_d = -1, math.inf
            for q in range(n):
                if q == p:
                    continue
                angle = math.atan2(
                    points[q, 1] - points[p, 1],
                    points[q, 0] - points[p, 0],
                ) % (2.0 * math.pi)
                if lo <= angle < hi and dist[p, q] < best_d:
                    best_d, best_q = dist[p, q], q
            if best_q >= 0:
                edges.add((min(p, best_q), max(p, best_q)))
    if not edges:
        return np.empty((0, 2), dtype=np.int64)
    return np.array(sorted(edges), dtype=np.int64)


# ── DGF (clean single-pass) ───────────────────────────────────────────────────

def _check_all_pairs_violation(
    edges_without: list[tuple[int, int]],
    dist: np.ndarray,
    w: np.ndarray,
    n: int,
    t: float,
    max_workers: int | None = None,
) -> bool:
    """
    Return True iff any pair (r, s) with r < s violates the t-spanner property
    in the graph defined by edges_without:
        w[r,s] * sp[r,s] > t * dist[r,s]   OR   sp[r,s] == inf

    Uses ThreadPoolExecutor with fast-break via threading.Event.
    """
    if n < 2:
        return False

    E = (
        np.array(edges_without, dtype=np.int64)
        if edges_without
        else np.empty((0, 2), dtype=np.int64)
    )
    adj = _build_adjacency(E, n, dist)
    sp = _floyd_warshall(adj)

    pairs = [(r, s) for r in range(n) for s in range(r + 1, n)]
    workers = max_workers or 4
    chunk_size = max(1, len(pairs) // workers)

    stop = threading.Event()
    found = threading.Event()

    def check_chunk(chunk: list[tuple[int, int]]) -> None:
        for r, s in chunk:
            if stop.is_set():
                return
            if np.isinf(sp[r, s]) or w[r, s] * sp[r, s] > t * dist[r, s]:
                found.set()
                stop.set()
                return

    chunks = [pairs[i : i + chunk_size] for i in range(0, len(pairs), chunk_size)]
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(check_chunk, c) for c in chunks]
        for f in as_completed(futs):
            if found.is_set():
                break

    return found.is_set()


def dgf_one_pass(
    edges: list[tuple[int, int]],
    dist: np.ndarray,
    n: int,
    t: float,
) -> tuple[list[tuple[int, int]], int]:
    """
    Single-pass Descending Greedy Filter (DGF).

    Process edges in non-ascending distance order (longest first).
    For each edge: tentatively remove it; if any pair (r, s) lacks a t-path
    in the remaining graph, reinsert it. ONE pass only — no rounds.

    Returns:
        (remaining_edges, n_removed)
    """
    if len(edges) <= 1:
        return list(edges), 0

    w = np.ones((n, n), dtype=np.float64)
    sel = sorted(edges, key=lambda e: -dist[e[0], e[1]])
    current = list(sel)
    removed = 0

    for edge in tqdm(sel, desc="dgf", unit="edge", leave=False):
        canon = (min(edge[0], edge[1]), max(edge[0], edge[1]))
        without = [e for e in current if (min(e[0], e[1]), max(e[0], e[1])) != canon]
        if not _check_all_pairs_violation(without, dist, w, n, t):
            current = without  # edge is redundant — remove it
            removed += 1
        # else: edge is necessary — leave it in current

    return current, removed


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


# ── Six algorithm wrappers ────────────────────────────────────────────────────

ALGO_CONFIG: dict[str, Any] = {"progress": True}


def algo_greedy(points: np.ndarray, t: float) -> np.ndarray:
    """Algorithm 1: standard greedy t-spanner."""
    return greedy_run(points, t, config=ALGO_CONFIG)


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
    stage1 = greedy_run(points, math.sqrt(t), config=ALGO_CONFIG)
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
    return greedy_run(points, t, E_input=yao_t_edges, config=ALGO_CONFIG)


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


def _save_results(out: Path, algo_results: dict[str, Any]) -> None:
    """Write results.json atomically (temp file + rename) so partial runs are safe."""
    tmp = out / "results.json.tmp"
    tmp.write_text(json.dumps(_to_jsonable(algo_results), indent=2))
    tmp.rename(out / "results.json")


def main() -> None:
    BASE_OUT = Path(f"results_final_experiment_n={N}")
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
