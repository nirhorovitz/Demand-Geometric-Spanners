"""Baseline profile artifact writer for SP-48 benchmark sweep."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Any

from core.json_utils import safe_dump


def write_baseline_artifact(
    out_dir: Path,
    profile: dict[str, Any],
    *,
    timestamp: datetime | None = None,
) -> Path:
    """
    Write baseline_profile.json and baseline_profile.md to out_dir.

    Args:
        out_dir: Directory to write artifacts (created if needed).
        profile: Baseline profile dict from run_baseline_profile.
        timestamp: Optional timestamp for report header.

    Returns:
        out_dir path.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ts = timestamp or datetime.now()
    ts_str = ts.strftime("%Y-%m-%dT%H:%M:%S")

    json_path = out_dir / "baseline_profile.json"
    with open(json_path, "w", encoding="utf-8") as f:
        safe_dump(profile, f, indent=2)

    md_lines = [
        "# SP-48 Baseline Profile",
        "",
        f"**Generated:** {ts_str}",
        f"**Problem type:** {profile.get('problem_type', 't')}",
        f"**Stretch t:** {profile.get('t')}",
        f"**Seed:** {profile.get('seed')}",
        f"**Repeats:** {profile.get('repeats')}",
        f"**Scales:** {profile.get('scales', [])}",
        f"**Algorithms:** {profile.get('algorithms', [])}",
        f"**Verification policy:** {profile.get('verification_policy', 'strict')}",
        "",
        "## Per-scale results",
        "",
    ]

    for n_str, scale_data in profile.get("per_scale", {}).items():
        md_lines.append(f"### n = {n_str}")
        md_lines.append("")
        for algo_id, run_data in scale_data.get("runs", {}).items():
            md_lines.append(f"- **{algo_id}**")
            md_lines.append(f"  - runtime_ms (median): {run_data.get('runtime_ms_median')}")
            md_lines.append(f"  - runtime_ms (p95): {run_data.get('runtime_ms_p95')}")
            md_lines.append(f"  - edge_count: {run_data.get('edge_count')}")
            md_lines.append(f"  - absolute_weight: {run_data.get('absolute_weight')}")
            md_lines.append(f"  - max_stretch: {run_data.get('max_stretch')}")
            md_lines.append(f"  - status: {run_data.get('status')}")
            if run_data.get("errors"):
                md_lines.append(f"  - errors: {run_data['errors']}")
            md_lines.append("")
        md_lines.append("")

    md_path = out_dir / "baseline_profile.md"
    with open(md_path, "w", encoding="utf-8") as f:
        f.write("\n".join(md_lines))

    return out_dir
