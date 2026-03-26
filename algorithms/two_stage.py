"""Reusable two-stage composition helper for chaining algorithms."""

from __future__ import annotations

import math
from typing import Any, Callable, Optional, Tuple

import numpy as np


def compute_stage1_params(t: float, config: Optional[dict[str, Any]] = None) -> Tuple[float, bool]:
    """
    Compute canonical stage1 parameters from (t, config) for reuse decisions.

    Returns:
        (t1, use_default_weight): t1 for stage1 stretch, use_default_weight for weight.
    """
    cfg = config or {}
    t1 = cfg.get("stage1_t")
    if t1 is None:
        t1 = math.sqrt(t)
    else:
        t1 = float(t1)
    use_default_weight = cfg.get("stage1_use_default_weight", True)
    return (t1, use_default_weight)


def _is_valid_stage1_output(arr: Any) -> bool:
    """Check if array is valid precomputed stage1 output: (m, 2) integer-like."""
    if arr is None:
        return False
    if not isinstance(arr, np.ndarray):
        return False
    if arr.ndim != 2 or arr.shape[1] != 2:
        return False
    if not np.issubdtype(arr.dtype, np.integer):
        return False
    return True


def run_two_stage(
    stage1: Callable[..., np.ndarray],
    stage2: Callable[..., np.ndarray],
    points: np.ndarray,
    t: float,
    *,
    E_input: Optional[np.ndarray] = None,
    weight: Optional[np.ndarray] = None,
    config: Optional[dict[str, Any]] = None,
    rng_seed: Optional[int] = None,
    stage1_output: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Compose two algorithms: stage-1 output edges are passed to stage-2 as E_input.

    Protocol-compatible inputs: points, t, E_input, weight, config, rng_seed.

    Optional stage1_output: when valid (m,2) integer-like array, skip stage1 and use
    it as E_input to stage2. Invalid or absent => run stage1 normally.

    Defaults:
        - Stage 1: t1 = sqrt(t) unless config["stage1_t"] provided
        - Stage 1: weight = ones unless config["stage1_use_default_weight"] is False
        - Stage 2: original t and original weight
        - Stage 2: E_input = stage-1 output edges

    Returns:
        Stage-2 output edges, shape (k, 2), same format as protocol.
    """
    cfg = config or {}

    # Stage 1 config
    t1, use_default_weight = compute_stage1_params(t, config)
    weight1: Optional[np.ndarray]
    if use_default_weight:
        weight1 = None  # ones convention
    else:
        weight1 = weight

    # Use precomputed stage1 when valid; otherwise run stage1
    if _is_valid_stage1_output(stage1_output):
        E_stage1 = np.asarray(stage1_output, dtype=np.int64).copy()
    else:
        # Run stage 1
        if hasattr(stage1, "run"):
            E_stage1 = stage1.run(
                points, t1,
                E_input=E_input,
                weight=weight1,
                config=config,
                rng_seed=rng_seed,
            )
        else:
            E_stage1 = stage1(
                points, t1,
                E_input=E_input,
                weight=weight1,
                config=config,
                rng_seed=rng_seed,
            )

    # Stage 2: original t, original weight, E_input = stage-1 output
    if hasattr(stage2, "run"):
        return stage2.run(
            points, t,
            E_input=E_stage1,
            weight=weight,
            config=config,
            rng_seed=rng_seed,
        )
    return stage2(
        points, t,
        E_input=E_stage1,
        weight=weight,
        config=config,
        rng_seed=rng_seed,
    )
