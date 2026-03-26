"""Run result dataclass for algorithm runs with metrics."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional

import numpy as np

from core.metrics import compute_metrics
from core.timing import perf_counter


@dataclass
class RunResult:
    """
    Result of an algorithm run with edges, metrics, max_stretch, and meta.
    max_stretch may be None when verification was skipped (large n).
    """

    edges: np.ndarray
    metrics: dict
    max_stretch: Optional[float]
    meta: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        self.edges = np.asarray(self.edges, dtype=np.int64)
        if self.edges.ndim == 1 and self.edges.size > 0:
            self.edges = self.edges.reshape(-1, 2)


def run_with_metrics(
    run_fn,
    points: np.ndarray,
    t: float,
    *,
    E_input: Optional[np.ndarray] = None,
    weight: Optional[np.ndarray] = None,
    config: Optional[dict] = None,
    rng_seed: Optional[int] = None,
    problem_type: str = "t",
    meta: Optional[dict] = None,
    skip_full_stretch_above_n: Optional[int] = None,
) -> RunResult:
    """
    Execute algorithm run function and wrap result with metrics.

    Args:
        run_fn: Callable matching AlgorithmProtocol.run signature.
        points: Point coordinates.
        t: Stretch parameter.
        E_input, weight, config, rng_seed: Passed to run_fn.
        problem_type: "t" or "tw".
        meta: Optional dict merged into RunResult.meta.
        skip_full_stretch_above_n: If n > this, skip O(n^3) stretch verification.

    Returns:
        RunResult with edges, metrics, max_stretch, meta.
    """
    t0 = perf_counter()
    edges = run_fn(
        points,
        t,
        E_input=E_input,
        weight=weight,
        config=config,
        rng_seed=rng_seed,
    )
    runtime_ms = (perf_counter() - t0) * 1000.0

    metrics = compute_metrics(
        edges=edges,
        points=points,
        weight=weight,
        t=t,
        runtime_ms=runtime_ms,
        problem_type=problem_type,
        skip_full_stretch_above_n=skip_full_stretch_above_n,
    )

    max_stretch = metrics.get("max_stretch")
    result_meta = dict(meta) if meta else {}
    result_meta.setdefault("problem_type", problem_type)
    result_meta.setdefault("t", t)
    result_meta.setdefault("n", points.shape[0])
    if "verification_note" in metrics:
        result_meta["verification_note"] = metrics["verification_note"]

    return RunResult(
        edges=edges,
        metrics=metrics,
        max_stretch=max_stretch,
        meta=result_meta,
    )
