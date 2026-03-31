"""Metric computation for spanner algorithm runs."""

from __future__ import annotations

from typing import Optional

import numpy as np

from core.mst import mst_weight

# Default tolerance for max_stretch <= t comparison
DEFAULT_STRETCH_TOLERANCE = 1e-9

# When n > this, full exact stretch verification (O(n^3)) is skipped.
# Metrics include verification_note; status is unknown (False for safety).
DEFAULT_SKIP_FULL_STRETCH_ABOVE_N = 2000


def _euclidean_distances(points: np.ndarray) -> np.ndarray:
    """Compute pairwise Euclidean distances for points (n, 2). Returns n x n matrix."""
    n = points.shape[0]
    if n == 0:
        return np.empty((0, 0), dtype=np.float64)
    diff = points[:, np.newaxis, :] - points[np.newaxis, :, :]
    return np.sqrt(np.sum(diff**2, axis=2))


def _build_adjacency(
    edges: np.ndarray,
    n: int,
    dist_matrix: np.ndarray,
) -> np.ndarray:
    """
    Build n x n adjacency matrix from edge set.
    dist_matrix[i,j] gives weight for edge (i,j).
    Non-edges get inf.
    """
    adj = np.full((n, n), np.inf, dtype=np.float64)
    np.fill_diagonal(adj, 0.0)
    if edges.size == 0:
        return adj
    for i, j in edges:
        i, j = int(i), int(j)
        w = dist_matrix[i, j]
        adj[i, j] = adj[j, i] = w
    return adj


from scipy.sparse.csgraph import shortest_path

def _floyd_warshall(adj: np.ndarray) -> np.ndarray:
    """All-pairs shortest paths via Floyd-Warshall (scipy). Returns n x n distance matrix."""
    # use scipy for C-level performance instead of Python triple-loops
    return shortest_path(adj, directed=False)


def apsp_shortest_path_after_removal(
    edges: np.ndarray,
    i: int,
    j: int,
    n: int,
    dist_matrix: np.ndarray,
) -> float:
    """
    Shortest path from i to j in graph G \\ {(i,j)} (graph with edge (i,j) removed).

    Uses Dijkstra from i for O(n^2) per call vs O(n^3) full Floyd-Warshall.
    Returns inf if i and j are disconnected after removal.

    Args:
        edges: (m, 2) edge array for current graph.
        i, j: Endpoints of the removed edge (0 <= i, j < n).
        n: Number of vertices.
        dist_matrix: n x n matrix of edge weights (lengths).

    Returns:
        Shortest path distance from i to j in G \\ {(i,j)}.
    """
    i, j = int(i), int(j)
    if i == j:
        return 0.0

    adj = _build_adjacency(edges, n, dist_matrix)
    # Remove edge (i,j) from adjacency
    adj[i, j] = np.inf
    adj[j, i] = np.inf

    dist = np.full(n, np.inf, dtype=np.float64)
    dist[i] = 0.0
    visited = np.zeros(n, dtype=bool)

    for _ in range(n):
        u = -1
        best = np.inf
        for v in range(n):
            if not visited[v] and dist[v] < best:
                best = dist[v]
                u = v
        if u < 0 or np.isinf(best):
            break
        if u == j:
            return dist[j]
        visited[u] = True
        for v in range(n):
            if not visited[v] and adj[u, v] < np.inf:
                cand = dist[u] + adj[u, v]
                if cand < dist[v]:
                    dist[v] = cand

    return dist[j] if not np.isinf(dist[j]) else np.inf


def apsp_add_edge(dist: np.ndarray, i: int, j: int, w: float) -> None:
    """
    Update APSP distance matrix in-place after adding one undirected weighted edge (i,j).

    For each pair (u,v), the new shortest path may use the new edge:
    dist[u,v] = min(dist[u,v], dist[u,i]+w+dist[j,v], dist[u,j]+w+dist[i,v]).

    Args:
        dist: n x n distance matrix (APSP for current graph), modified in-place.
        i, j: Endpoints of the new edge (0 <= i, j < n).
        w: Weight of the new edge (must be finite and non-negative).
    """
    n = dist.shape[0]
    i, j = int(i), int(j)
    if i == j:
        return
    for u in range(n):
        for v in range(n):
            cand1 = dist[u, i] + w + dist[j, v]
            cand2 = dist[u, j] + w + dist[i, v]
            new_val = min(dist[u, v], cand1, cand2)
            dist[u, v] = new_val


def _compute_max_stretch(
    edges: np.ndarray,
    n: int,
    dist_matrix: np.ndarray,
) -> float:
    """
    Compute max stretch over all pairs in the same connected component.
    stretch(u,v) = shortest_path(u,v) / direct_distance(u,v).
    Returns inf if disconnected or zero-edge (no path between some pair).
    Returns 0.0 only when n <= 1.
    """
    if n <= 1:
        return 0.0
    if edges.size == 0:
        return np.inf  # Disconnected: no path between any pair

    adj = _build_adjacency(edges, n, dist_matrix)
    shortest = _floyd_warshall(adj)

    max_stretch = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            direct = dist_matrix[i, j]
            path_len = shortest[i, j]
            if np.isinf(path_len):
                return np.inf  # Disconnected
            if direct > 0:
                stretch = path_len / direct
                if stretch > max_stretch:
                    max_stretch = stretch
    return max_stretch


def _stretch_verification_skipped_note(n: int, threshold: int) -> str:
    """Build metadata note when full stretch verification is skipped."""
    return (
        f"Full exact stretch verification skipped for n={n} > {threshold}. "
        "Verification is O(n^3); use smaller n for exact metrics."
    )


def _compute_degrees(edges: np.ndarray, n: int) -> tuple[int, float]:
    """Return (max_degree, average_degree). For zero edges: (0, 0.0)."""
    if n == 0 or edges.size == 0:
        return 0, 0.0
    deg = np.zeros(n, dtype=np.int64)
    for i, j in edges:
        deg[int(i)] += 1
        deg[int(j)] += 1
    return int(deg.max()), float(deg.mean())


def _absolute_weight(edges: np.ndarray, dist_matrix: np.ndarray) -> float:
    """Sum of edge weights. Zero edges => 0.0."""
    if edges.size == 0:
        return 0.0
    total = 0.0
    for i, j in edges:
        total += dist_matrix[int(i), int(j)]
    return total


def compute_metrics(
    edges: np.ndarray,
    points: np.ndarray,
    weight: Optional[np.ndarray],
    t: float,
    runtime_ms: float,
    problem_type: str = "t",
    stretch_tolerance: float = DEFAULT_STRETCH_TOLERANCE,
    skip_full_stretch_above_n: Optional[int] = None,
) -> dict:
    """
    Compute required metrics for an algorithm run.

    Args:
        edges: Output edge set, shape (k, 2).
        points: Point coordinates, shape (n, 2).
        weight: n x n weight matrix for tw mode; None for t mode.
        t: Stretch threshold; status valid iff max_stretch <= t (within tolerance).
        runtime_ms: Elapsed runtime in milliseconds.
        problem_type: "t" (Euclidean) or "tw" (weighted).
        stretch_tolerance: Tolerance for max_stretch <= t comparison.
        skip_full_stretch_above_n: If n > this, skip O(n^3) stretch verification.
            Uses DEFAULT_SKIP_FULL_STRETCH_ABOVE_N when None.

    Returns:
        Dict with: runtime_ms, edge_count, status, absolute_weight,
        relative_weight_to_mst, max_degree, average_degree, max_stretch.
        When stretch verification is skipped: verification_note, max_stretch=None.
    """
    print("        [compute_metrics] Starting...", flush=True)
    n = points.shape[0]
    edges = np.asarray(edges, dtype=np.int64)
    if edges.ndim == 1 and edges.size > 0:
        edges = edges.reshape(-1, 2)

    threshold = (
        skip_full_stretch_above_n
        if skip_full_stretch_above_n is not None
        else DEFAULT_SKIP_FULL_STRETCH_ABOVE_N
    )
    skip_stretch = n > threshold

    print(f"        [compute_metrics] n={n}, threshold={threshold}, skip_stretch={skip_stretch}", flush=True)

    if problem_type == "t":
        dist_matrix = _euclidean_distances(points)
    else:
        # tw mode: use weight matrix as distances
        if weight is None:
            weight = np.ones((n, n), dtype=np.float64)
        dist_matrix = np.asarray(weight, dtype=np.float64)

    if skip_stretch:
        print("        [compute_metrics] Skipping stretch verification...", flush=True)
        max_stretch = None
        verification_note = _stretch_verification_skipped_note(n, threshold)
        status = False  # Unknown; conservative
    else:
        print("        [compute_metrics] Computing max stretch...", flush=True)
        max_stretch = _compute_max_stretch(edges, n, dist_matrix)
        print("        [compute_metrics] Max stretch computed.", flush=True)
        verification_note = None

    print("        [compute_metrics] Computing absolute weight...", flush=True)
    abs_w = _absolute_weight(edges, dist_matrix)
    print("        [compute_metrics] Computing MST weight...", flush=True)
    mst_w = mst_weight(points, weight, problem_type)
    print("        [compute_metrics] Computing degrees...", flush=True)
    max_deg, avg_deg = _compute_degrees(edges, n)

    print("        [compute_metrics] Evaluating final metrics...", flush=True)

    # relative_weight_to_mst: abs_w / mst_w. Handle zero MST (n<=1) and zero edges.
    if mst_w <= 0:
        rel_to_mst = 0.0 if abs_w <= 0 else np.inf
    else:
        rel_to_mst = abs_w / mst_w
    print("        [compute_metrics] Relative weight to MST computed.", flush=True)

    # status: valid iff max_stretch <= t (with tolerance). When skipped, False.
    if not skip_stretch:
        if np.isinf(max_stretch):
            status = False
        else:
            status = max_stretch <= t or np.isclose(
                max_stretch, t, rtol=0, atol=stretch_tolerance
            )
        status = bool(status)
    print("        [compute_metrics] Status computed.", flush=True)
    out: dict = {
        "runtime_ms": runtime_ms,
        "edge_count": int(edges.shape[0]),
        "status": status,
        "absolute_weight": abs_w,
        "relative_weight_to_mst": rel_to_mst,
        "max_degree": max_deg,
        "average_degree": avg_deg,
        "max_stretch": (
            (max_stretch if not np.isinf(max_stretch) else float("inf"))
            if max_stretch is not None
            else None
        ),
    }
    if verification_note is not None:
        out["verification_note"] = verification_note
    print("        [compute_metrics] Returning results.", flush=True)
    return out
