"""Shared algorithm protocol and base utilities."""

import sys
from typing import Any, Iterable, Iterator, Optional, Protocol, runtime_checkable, TypeVar

import numpy as np

from core.graph_utils import get_candidate_edges

try:
    from tqdm.auto import tqdm
except ImportError:  # pragma: no cover - dependency fallback
    tqdm = None  # type: ignore[assignment]

T = TypeVar("T")


@runtime_checkable
class AlgorithmProtocol(Protocol):
    """Protocol for spanner algorithms."""

    def run(
        self,
        points: np.ndarray,
        t: float,
        *,
        E_input: Optional[np.ndarray] = None,
        weight: Optional[np.ndarray] = None,
        config: Optional[dict[str, Any]] = None,
        rng_seed: Optional[int] = None,
    ) -> np.ndarray:
        """
        Run the algorithm on a point set.

        Args:
            points: Point coordinates, shape (n, 2).
            t: Spanner stretch parameter.
            E_input: Candidate edge set, shape (m, 2). None => full graph.
            weight: n x n weight matrix for tw mode. None => ones convention.
            config: Optional algorithm-specific config.
                Common progress keys:
                - "progress": bool (force enable/disable tqdm bars)
                - "progress_min_total": int (min iterations before showing)
            rng_seed: Optional random seed for reproducibility.

        Returns:
            Selected edges, shape (k, 2), each row (i, j).
        """
        ...


def resolve_weight(weight: Optional[np.ndarray], n: int) -> np.ndarray:
    """
    Resolve weight matrix. weight=None => ones convention.

    Returns:
        n x n weight matrix.
    """
    if weight is None:
        return np.ones((n, n), dtype=np.float64)
    return np.asarray(weight, dtype=np.float64)


def resolve_candidates(points: np.ndarray, E_input: Optional[np.ndarray]) -> np.ndarray:
    """
    Resolve candidate edge set for algorithm.

    E_input=None => full graph candidates.
    E_input provided => validated edge set for chaining.

    Returns:
        Edge array of shape (m, 2).
    """
    n = points.shape[0]
    return get_candidate_edges(n, E_input)


APSP_MODE_FULL = "full"
APSP_MODE_INCREMENTAL = "incremental"


def resolve_apsp_mode(config: Optional[dict[str, Any]]) -> str:
    """
    Resolve apsp_mode from config for incremental vs full Floyd-Warshall.

    config["apsp_mode"] in {"full", "incremental"}.
    Default: "incremental" for runtime optimization; "full" for fallback/validation.
    """
    if config is None:
        return APSP_MODE_INCREMENTAL
    mode = config.get("apsp_mode", APSP_MODE_INCREMENTAL)
    if mode in (APSP_MODE_FULL, APSP_MODE_INCREMENTAL):
        return mode
    return APSP_MODE_INCREMENTAL


def progress_iter(
    iterable: Iterable[T],
    *,
    total: Optional[int] = None,
    desc: str = "",
    config: Optional[dict[str, Any]] = None,
) -> Iterator[T]:
    """
    Wrap iterable with tqdm based on protocol config.

    Defaults to showing progress in interactive terminals and hiding it in
    non-interactive contexts (e.g., test runs, CI logs).
    """
    cfg = config or {}
    enabled = cfg.get("progress")
    if enabled is None:
        enabled = sys.stderr.isatty()
    min_total = int(cfg.get("progress_min_total", 100))

    if (
        bool(enabled)
        and tqdm is not None
        and (total is None or total >= min_total)
    ):
        return tqdm(iterable, total=total, desc=desc, leave=False)
    return iter(iterable)
