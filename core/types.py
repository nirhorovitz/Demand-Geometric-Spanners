"""Core protocol dataclasses and types for Spanner experiments."""

from dataclasses import dataclass
from typing import Literal, Optional

import numpy as np


# --- Point set protocol ---

SourceType = Literal["random", "fixed"]


@dataclass
class PointSet:
    """Point set matching the protocol: set_id, points, n, source_type, generator_meta."""

    set_id: str
    points: np.ndarray  # shape (n, 2), each row is (x, y)
    n: int
    source_type: SourceType
    generator_meta: Optional[dict] = None

    def __post_init__(self):
        self.points = np.asarray(self.points, dtype=np.float64)
        if self.generator_meta is None:
            self.generator_meta = {}


# --- Experiment config ---

ProblemType = Literal["t", "tw"]


@dataclass
class ExperimentConfig:
    """Experiment configuration: problem_type, t, n, algorithms, seed."""

    problem_type: ProblemType
    t: float
    n: int
    algorithms: list[str]
    seed: Optional[int] = None
    weight: Optional[np.ndarray] = None  # for tw: n x n matrix; None = ones convention
