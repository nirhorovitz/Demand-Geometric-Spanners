"""
Filtered greedy: greedy stage1 (defaults) -> filter stage2 (original t/weight).

Stage 1: Greedy t-spanner with helper defaults (t1=sqrt(t), weight=ones).
Stage 2: Filter remove/reinsert validation with original t and weight.
"""

from __future__ import annotations

from typing import Any, Optional

from algorithms import filter_algorithm
from algorithms import greedy
from algorithms.registry import register
from algorithms.two_stage import run_two_stage


@register("filtered_greedy")
def run(
    points: np.ndarray,
    t: float,
    *,
    E_input: Optional[np.ndarray] = None,
    weight: Optional[np.ndarray] = None,
    config: Optional[dict[str, Any]] = None,
    rng_seed: Optional[int] = None,
) -> np.ndarray:
    """
    Filtered greedy: greedy(defaults) -> filter(original t/weight).

    Stage 1: greedy with t1=sqrt(t), weight=ones.
    Stage 2: filter with original t and weight on stage-1 output edges.
    """
    precomputed = config.get("_stage1_output") if config else None
    return run_two_stage(
        greedy.run,
        filter_algorithm.run,
        points,
        t,
        E_input=E_input,
        weight=weight,
        config=config,
        rng_seed=rng_seed,
        stage1_output=precomputed,
    )
