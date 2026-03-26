"""JSON serialization utilities for experiment outputs."""

from __future__ import annotations

import json
from typing import Any

import numpy as np


def _json_default(obj: Any) -> Any:
    """Convert non-JSON-serializable objects for json.dumps."""
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.integer, np.int64, np.int32)):
        return int(obj)
    if isinstance(obj, (np.floating, np.float64, np.float32)):
        f = float(obj)
        if np.isinf(f):
            return "Infinity" if f > 0 else "-Infinity"
        if np.isnan(f):
            return "NaN"
        return f
    if isinstance(obj, float):
        if np.isinf(obj):
            return "Infinity" if obj > 0 else "-Infinity"
        if np.isnan(obj):
            return "NaN"
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")


def _to_json_safe(obj: Any) -> Any:
    """Recursively convert to JSON-serializable form (handles float inf/nan, numpy)."""
    if isinstance(obj, np.ndarray):
        return _to_json_safe(obj.tolist())
    if isinstance(obj, (np.integer, np.int64, np.int32)):
        return int(obj)
    if isinstance(obj, (np.floating, np.float64, np.float32)):
        f = float(obj)
        if np.isinf(f):
            return "Infinity" if f > 0 else "-Infinity"
        if np.isnan(f):
            return "NaN"
        return f
    if isinstance(obj, float):
        if np.isinf(obj):
            return "Infinity" if obj > 0 else "-Infinity"
        if np.isnan(obj):
            return "NaN"
    if isinstance(obj, dict):
        return {k: _to_json_safe(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_to_json_safe(v) for v in obj]
    return obj


def sanitize_metadata(meta: dict[str, Any]) -> dict[str, Any]:
    """
    Convert metadata to JSON-serializable form.

    Replaces non-serializable values (numpy, Path, etc.) with string repr
    or drops them. Never raises.
    """
    out: dict[str, Any] = {}
    for k, v in meta.items():
        try:
            json.dumps(v, default=_json_default)
            out[k] = v
        except (TypeError, ValueError):
            try:
                out[k] = str(v)
            except Exception:
                out[k] = "<non-serializable>"
    return out


def safe_dumps(obj: Any, **kwargs: Any) -> str:
    """
    Serialize to JSON with numpy and common types supported.

    Uses _json_default for numpy arrays; _to_json_safe for float inf/nan in nested structures.
    """
    preprocessed = _to_json_safe(obj)
    return json.dumps(preprocessed, default=_json_default, **kwargs)


def safe_dump(obj: Any, fp: Any, **kwargs: Any) -> None:
    """Write JSON to file with numpy support."""
    preprocessed = _to_json_safe(obj)
    json.dump(preprocessed, fp, default=_json_default, **kwargs)
