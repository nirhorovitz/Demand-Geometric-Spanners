"""Schema validation hook for manifest and results before writing."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from jsonschema import Draft7Validator, ValidationError


def _load_schema(name: str) -> dict[str, Any]:
    """Load JSON schema from schemas/ directory."""
    schema_dir = Path(__file__).resolve().parents[1] / "schemas"
    path = schema_dir / f"{name}.schema.json"
    import json
    with open(path, encoding="utf-8") as f:
        return json.load(f)


_MANIFEST_SCHEMA = None
_RESULTS_SCHEMA = None


def _get_manifest_schema() -> dict[str, Any]:
    global _MANIFEST_SCHEMA
    if _MANIFEST_SCHEMA is None:
        _MANIFEST_SCHEMA = _load_schema("manifest")
    return _MANIFEST_SCHEMA


def _get_results_schema() -> dict[str, Any]:
    global _RESULTS_SCHEMA
    if _RESULTS_SCHEMA is None:
        _RESULTS_SCHEMA = _load_schema("results")
    return _RESULTS_SCHEMA


def validate_manifest(payload: dict[str, Any]) -> None:
    """
    Validate manifest payload against manifest schema.

    Raises:
        jsonschema.ValidationError: With clear message if invalid.
    """
    validator = Draft7Validator(_get_manifest_schema())
    errors = list(validator.iter_errors(payload))
    if errors:
        msg = "; ".join(
            f"{e.json_path or 'root'}: {e.message}" for e in errors
        )
        raise ValidationError(msg)


def validate_results(payload: dict[str, Any]) -> None:
    """
    Validate results payload against results schema.

    Raises:
        jsonschema.ValidationError: With clear message if invalid.
    """
    validator = Draft7Validator(_get_results_schema())
    errors = list(validator.iter_errors(payload))
    if errors:
        msg = "; ".join(
            f"{e.json_path or 'root'}: {e.message}" for e in errors
        )
        raise ValidationError(msg)
