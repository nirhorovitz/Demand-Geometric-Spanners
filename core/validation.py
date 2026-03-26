"""Validation for point sets and weight matrices."""

from typing import Optional

import numpy as np

from core.types import PointSet


def validate_point_set(ps: PointSet) -> None:
    """Validate point set: n consistency, points in [0,1]x[0,1], protocol compliance."""
    if ps.points.shape[0] != ps.n:
        raise ValueError(
            f"Point set n mismatch: declared n={ps.n}, actual points count={ps.points.shape[0]}"
        )
    if ps.points.ndim != 2 or ps.points.shape[1] != 2:
        raise ValueError(
            f"Points must be 2D array of shape (n, 2), got shape {ps.points.shape}"
        )
    if np.any(ps.points < 0) or np.any(ps.points > 1):
        raise ValueError(
            "All point coordinates must be in [0, 1]"
        )


def validate_weight_matrix(
    weight: Optional[np.ndarray],
    n: int,
) -> np.ndarray:
    """
    Validate weight matrix for tw problem.
    - shape n x n
    - symmetric
    - values in [0, 1]
    - weight=None resolves to all-ones matrix convention.

    Returns the validated (or default) weight matrix.
    """
    if weight is None:
        return np.ones((n, n), dtype=np.float64)

    w = np.asarray(weight, dtype=np.float64)
    if w.shape != (n, n):
        raise ValueError(
            f"Weight matrix must be n x n (n={n}), got shape {w.shape}"
        )
    if not np.allclose(w, w.T):
        raise ValueError("Weight matrix must be symmetric")
    if np.any(w < 0) or np.any(w > 1):
        raise ValueError("Weight matrix values must be in [0, 1]")
    return w
