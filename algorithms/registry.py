"""Algorithm registry: dispatch by id, clear errors on unknown."""

from typing import Any, Callable, Dict, Optional, Union

import numpy as np

from algorithms.base import AlgorithmProtocol
from core.run_result import RunResult, run_with_metrics


# Callable-compatible protocol signature (supports keyword-only protocol args)
RunSignature = Callable[..., np.ndarray]


_REGISTRY: Dict[str, Union[AlgorithmProtocol, RunSignature]] = {}


def register(algorithm_id: str) -> Callable:
    """Decorator to register an algorithm by id."""

    def _decorator(fn: Union[AlgorithmProtocol, RunSignature]):
        _REGISTRY[algorithm_id] = fn
        return fn

    return _decorator


def get(algorithm_id: str) -> Union[AlgorithmProtocol, RunSignature]:
    """
    Get algorithm by id.

    Raises:
        KeyError: With clear message if algorithm_id is unknown.
    """
    if algorithm_id not in _REGISTRY:
        available = ", ".join(list_algorithms()) or "(none)"
        raise KeyError(
            f"Unknown algorithm: {algorithm_id!r}. "
            f"Available: {available}"
        )
    return _REGISTRY[algorithm_id]


def run_algorithm(
    algorithm_id: str,
    points: np.ndarray,
    t: float,
    *,
    E_input: Optional[np.ndarray] = None,
    weight: Optional[np.ndarray] = None,
    config: Optional[dict] = None,
    rng_seed: Optional[int] = None,
) -> np.ndarray:
    """
    Run algorithm by id with full protocol signature.

    Dispatches to the registered algorithm's run method.
    """
    algo = get(algorithm_id)
    if hasattr(algo, "run"):
        return algo.run(
            points, t,
            E_input=E_input,
            weight=weight,
            config=config,
            rng_seed=rng_seed,
        )
    # Plain callable: invoke with full protocol signature
    return algo(
        points, t,
        E_input=E_input,
        weight=weight,
        config=config,
        rng_seed=rng_seed,
    )


# Canonical order for algorithm listing. Unknown IDs are appended sorted after this block.
_CANONICAL_ORDER = [
    "greedy",
    "filtered_greedy",
    "filter",
    "repaired_greedy",
    "repair",
    "assignment_greedy",
    "double_greedy",
    "bucket_greedy",
    "weighted_ordered_greedy",
]


def list_algorithms() -> list[str]:
    """Return algorithms in canonical order; unknown IDs appended sorted."""
    known = [aid for aid in _CANONICAL_ORDER if aid in _REGISTRY]
    unknown = sorted(k for k in _REGISTRY.keys() if k not in _CANONICAL_ORDER)
    return known + unknown


def run_algorithm_with_metrics(
    algorithm_id: str,
    points: np.ndarray,
    t: float,
    *,
    E_input: Optional[np.ndarray] = None,
    weight: Optional[np.ndarray] = None,
    config: Optional[dict] = None,
    rng_seed: Optional[int] = None,
    problem_type: str = "t",
    meta: Optional[dict[str, Any]] = None,
    skip_full_stretch_above_n: Optional[int] = None,
) -> RunResult:
    """
    Run algorithm by id and return RunResult with edges, metrics, max_stretch, meta.
    """
    algo = get(algorithm_id)

    def _run(points: np.ndarray, t: float, **kwargs) -> np.ndarray:
        if hasattr(algo, "run"):
            return algo.run(points, t, **kwargs)
        return algo(points, t, **kwargs)

    result_meta = dict(meta) if meta else {}
    result_meta.setdefault("algorithm_id", algorithm_id)

    return run_with_metrics(
        _run,
        points,
        t,
        E_input=E_input,
        weight=weight,
        config=config,
        rng_seed=rng_seed,
        problem_type=problem_type,
        meta=result_meta,
        skip_full_stretch_above_n=skip_full_stretch_above_n,
    )
