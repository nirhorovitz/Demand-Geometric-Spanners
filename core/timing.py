"""Runtime instrumentation for algorithm and phase timing."""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Any, Callable, TypeVar

T = TypeVar("T")


@dataclass
class TimingResult:
    """
    Result of a timed operation with optional phase breakdown.

    Attributes:
        runtime_ms: Total elapsed time in milliseconds.
        phases: Optional dict of phase name -> ms for breakdown.
        meta: Optional metadata (e.g., verification_skipped note).
    """

    runtime_ms: float
    phases: dict[str, float] = field(default_factory=dict)
    meta: dict[str, Any] = field(default_factory=dict)

    def to_metrics_metadata(self) -> dict[str, Any]:
        """Merge timing metadata for inclusion in run metrics."""
        out: dict[str, Any] = {"runtime_ms": self.runtime_ms}
        if self.phases:
            out["phase_timings_ms"] = self.phases
        out.update(self.meta)
        return out


def measure_fn(fn: Callable[..., T], *args: Any, **kwargs: Any) -> tuple[T, float]:
    """
    Run a callable and return (result, runtime_ms).

    Args:
        fn: Callable to execute.
        *args, **kwargs: Passed to fn.

    Returns:
        Tuple of (fn(*args, **kwargs), runtime_ms).
    """
    t0 = time.perf_counter()
    result = fn(*args, **kwargs)
    runtime_ms = (time.perf_counter() - t0) * 1000.0
    return result, runtime_ms


def perf_counter() -> float:
    """
    Return current time from perf_counter for manual timing.

    Usage:
        t0 = perf_counter()
        # ... work ...
        runtime_ms = (perf_counter() - t0) * 1000.0
    """
    return time.perf_counter()
