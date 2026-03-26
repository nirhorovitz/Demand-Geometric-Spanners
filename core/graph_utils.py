"""Graph utilities for algorithm chaining and edge set handling."""

from typing import Optional

import numpy as np


def full_graph_candidates(n: int) -> np.ndarray:
    """
    Return all undirected edges as candidate set for n nodes.

    Format: array of shape (m, 2), each row (i, j) with i < j.
    No self-loops, no duplicates.
    """
    rows = []
    for i in range(n):
        for j in range(i + 1, n):
            rows.append([i, j])
    return np.array(rows, dtype=np.int64) if rows else np.empty((0, 2), dtype=np.int64)


def validate_edge_set(E: np.ndarray, n: int) -> np.ndarray:
    """
    Validate and normalize edge set for n nodes.

    - Raises ValueError if any node index is out of range [0, n).
    - Removes self-loops and duplicate edges.
    - Returns normalized edges as (m, 2) with i < j per row.

    Returns:
        Cleaned edge array of shape (m, 2).
    """
    E = np.asarray(E, dtype=np.int64)
    if E.size == 0:
        return np.empty((0, 2), dtype=np.int64)
    if E.ndim != 2 or E.shape[1] != 2:
        raise ValueError(f"Edge set must be (m, 2), got shape {E.shape}")

    # Out-of-range check
    if np.any(E < 0) or np.any(E >= n):
        raise ValueError(
            f"Edge indices must be in [0, n) for n={n}; "
            f"got min={E.min()}, max={E.max()}"
        )

    # Normalize: (i, j) -> (min, max) to canonicalize
    E = np.sort(E, axis=1)
    # Drop self-loops
    E = E[E[:, 0] != E[:, 1]]
    # Drop duplicates
    E = np.unique(E, axis=0)
    return E


def get_candidate_edges(
    n: int,
    E_input: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Resolve candidate edge set for algorithm input.

    - E_input=None => full graph candidates.
    - E_input provided => validated and cleaned edge set (for chaining).

    Returns:
        Edge array of shape (m, 2).
    """
    if E_input is None:
        return full_graph_candidates(n)
    return validate_edge_set(E_input, n)
