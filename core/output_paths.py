"""Canonical output path builder for experiment results."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Literal

ProblemType = Literal["t", "tw"]


def build_results_dir_name(
    problem_type: ProblemType,
    set_id: str,
    n: int,
    t: float,
    algorithms: list[str],
    *,
    timestamp: datetime | None = None,
) -> str:
    """
    Build canonical output directory name.

    Format: problem-{t|tw}__{set_label}__t={t}__algos={k}__ts={MMDD-HHMMSS}

    Args:
        problem_type: "t" or "tw".
        set_id: Point set identifier.
        n: Number of points.
        t: Stretch parameter.
        algorithms: List of algorithm ids (k = len).
        timestamp: Optional datetime; defaults to now.

    Returns:
        Directory name string (no path separators).
    """
    ts = timestamp or datetime.now()
    ts_str = ts.strftime("%m%d-%H%M%S")
    k = len(algorithms)

    raw_set = str(set_id)
    lower_set = raw_set.lower()

    # Custom point-set label formatting requested by user:
    # random -> n=## ; fixed -> fixed_set_of_points:id=##
    if lower_set.startswith("random"):
        set_label = f"n={n}"
    elif lower_set.startswith("fixed"):
        fixed_id = raw_set[len("fixed"):].lstrip("-_:")
        if not fixed_id:
            fixed_id = raw_set
        set_label = f"fixed_set_of_points:id={fixed_id}"
    else:
        safe_set = raw_set.replace("/", "_").replace("\\", "_")
        set_label = safe_set

    return f"problem-{problem_type}__{set_label}__t={t}__algos={k}__ts={ts_str}"


def resolve_output_path(
    base_dir: Path,
    problem_type: ProblemType,
    set_id: str,
    n: int,
    t: float,
    algorithms: list[str],
    *,
    timestamp: datetime | None = None,
) -> Path:
    """
    Resolve output directory path, handling collisions.

    If the canonical path already exists (e.g., same second), appends
    suffix _1, _2, ... until a non-existent path is found.

    Args:
        base_dir: Root results directory.
        problem_type, set_id, n, t, algorithms: Same as build_results_dir_name.
        timestamp: Optional datetime.

    Returns:
        Path to a non-existent directory ready for creation.
    """
    base_dir = Path(base_dir)
    base_name = build_results_dir_name(
        problem_type, set_id, n, t, algorithms, timestamp=timestamp
    )
    path = base_dir / base_name
    if not path.exists():
        return path

    # Collision: append _1, _2, ...
    idx = 1
    while True:
        candidate = base_dir / f"{base_name}_{idx}"
        if not candidate.exists():
            return candidate
        idx += 1
