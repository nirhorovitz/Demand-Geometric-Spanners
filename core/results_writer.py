"""Results writer for experiment outputs with schema validation."""

from __future__ import annotations

import platform
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Literal

import numpy as np

from core.json_utils import safe_dump, sanitize_metadata
from core.output_paths import resolve_output_path
from core.schema_validation import validate_manifest, validate_results


def run_result_to_dict(result: Any) -> dict[str, Any]:
    """
    Convert RunResult (or similar) to dict for write_experiment_results.

    Handles objects with attributes: edges, metrics, max_stretch, meta.
    """
    meta = getattr(result, "meta", None) or {}
    return {
        "algorithm_id": meta.get("algorithm_id", ""),
        "edges": getattr(result, "edges", None),
        "metrics": getattr(result, "metrics", None),
        "max_stretch": getattr(result, "max_stretch", None),
        "meta": meta,
    }


ProblemType = Literal["t", "tw"]


def _gather_runtime_env() -> dict[str, Any]:
    """Gather runtime environment metadata (JSON-serializable)."""
    return sanitize_metadata({
        "python_version": sys.version,
        "platform": platform.platform(),
        "python_implementation": platform.python_implementation(),
    })


def _run_to_dict(
    algorithm_id: str,
    edges: np.ndarray | None,
    metrics: dict | None,
    max_stretch: float | None,
    meta: dict[str, Any] | None,
    error: str | None = None,
) -> dict[str, Any]:
    """Convert a single algorithm run (success or failure) to JSON-serializable dict."""
    out: dict[str, Any] = {"algorithm_id": algorithm_id}
    if error is not None:
        out["error"] = error
        out["edges"] = []
        out["metrics"] = {}
        out["max_stretch"] = float("inf") if max_stretch is None else max_stretch
        out["meta"] = sanitize_metadata(meta or {})
        return out

    edges_arr = np.asarray(edges, dtype=np.int64)
    if edges_arr.ndim == 1 and edges_arr.size > 0:
        edges_arr = edges_arr.reshape(-1, 2)
    out["edges"] = edges_arr.tolist()
    out["metrics"] = metrics or {}
    # max_stretch may be None when verification was skipped (large n)
    out["max_stretch"] = max_stretch
    out["meta"] = sanitize_metadata(meta or {})
    return out


def write_experiment_results(
    *,
    base_dir: Path,
    problem_type: ProblemType,
    set_id: str,
    n: int,
    t: float,
    algorithms: list[str],
    points: np.ndarray,
    run_results: list[dict[str, Any]],
    seed: int | None = None,
    timestamp: datetime | None = None,
) -> Path:
    """
    Write manifest.json and results.json to a canonical output directory.

    Args:
        base_dir: Root results directory.
        problem_type: "t" or "tw".
        set_id: Point set identifier.
        n: Number of points.
        t: Stretch parameter.
        algorithms: List of algorithm ids (order preserved).
        points: Point coordinates, shape (n, 2).
        run_results: List of per-algorithm results. Each element is either:
            - RunResult-like: {"algorithm_id", "edges", "metrics", "max_stretch", "meta"}
            - Failure: {"algorithm_id", "error", "meta"} (edges/metrics may be absent)
        seed: Optional random seed.
        timestamp: Optional datetime for path naming.

    Returns:
        Path to the created output directory.

    Raises:
        jsonschema.ValidationError: If payloads fail schema validation.
    """
    base_dir = Path(base_dir)
    out_path = resolve_output_path(
        base_dir, problem_type, set_id, n, t, algorithms, timestamp=timestamp
    )
    out_path.mkdir(parents=True, exist_ok=False)

    ts = timestamp or datetime.now()
    ts_str = ts.strftime("%Y-%m-%dT%H:%M:%S")
    experiment_id = out_path.name

    config = {
        "problem_type": problem_type,
        "set_id": set_id,
        "n": n,
        "t": t,
        "algorithms": algorithms,
        "seed": seed,
    }

    # Build manifest
    manifest = {
        "experiment_id": experiment_id,
        "timestamp": ts_str,
        "output_dir": str(out_path),
        "config": config,
        "runtime_env": _gather_runtime_env(),
    }
    validate_manifest(manifest)

    # Build results with all points and all edges per algorithm
    points_list = np.asarray(points, dtype=np.float64)
    if points_list.ndim == 1:
        points_list = points_list.reshape(-1, 2)
    points_list = points_list.tolist()

    algorithm_runs = []
    for r in run_results:
        algo_id = r.get("algorithm_id", "")
        err = r.get("error")
        algorithm_runs.append(
            _run_to_dict(
                algorithm_id=algo_id,
                edges=r.get("edges"),
                metrics=r.get("metrics"),
                max_stretch=r.get("max_stretch"),
                meta=r.get("meta"),
                error=err,
            )
        )

    results_payload = {
        "points": points_list,
        "algorithm_runs": algorithm_runs,
        "config": config,
    }
    validate_results(results_payload)

    # Write files
    with open(out_path / "manifest.json", "w", encoding="utf-8") as f:
        safe_dump(manifest, f, indent=2)

    with open(out_path / "results.json", "w", encoding="utf-8") as f:
        safe_dump(results_payload, f, indent=2)

    return out_path
