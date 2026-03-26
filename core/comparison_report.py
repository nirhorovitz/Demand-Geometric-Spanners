"""SP-59: Before/after benchmark comparison report for filter, repair, weighted_greedy_spanner trio."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Any

from compressions.common import SP48_TRIO_ALGORITHMS
from core.json_utils import safe_dump


# Canonical trio for SP-48/SP-59. Filter to these only in reports.
TRIO_ALGORITHMS = tuple(SP48_TRIO_ALGORITHMS)


def filter_trio_only(per_scale: dict[str, Any]) -> dict[str, Any]:
    """
    Filter per_scale runs to only include filter, repair, greedy (trio).

    Args:
        per_scale: Dict mapping scale str -> {"n", "set_id", "runs": {algo_id: run_data}}.

    Returns:
        New dict with same structure but runs restricted to trio algorithms.
    """
    trio_set = set(TRIO_ALGORITHMS)
    out: dict[str, Any] = {}
    for n_str, scale_data in per_scale.items():
        runs = scale_data.get("runs", {})
        filtered_runs = {aid: data for aid, data in runs.items() if aid in trio_set}
        out[n_str] = {
            **scale_data,
            "runs": filtered_runs,
        }
    return out


def _extract_run_metrics(run_data: dict[str, Any]) -> dict[str, Any]:
    """Extract runtime_ms_median, max_stretch, status, errors from run_data."""
    return {
        "runtime_ms_median": run_data.get("runtime_ms_median"),
        "max_stretch": run_data.get("max_stretch"),
        "status": run_data.get("status"),
        "errors": run_data.get("errors"),
    }


def _compute_deltas(
    before: dict[str, Any],
    after: dict[str, Any],
) -> dict[str, Any]:
    """
    Compute runtime deltas (absolute and %) between before and after.

    Returns dict with runtime_delta_ms, runtime_delta_pct, or None if either failed.
    """
    rb = before.get("runtime_ms_median")
    ra = after.get("runtime_ms_median")
    if rb is None or ra is None:
        return {"runtime_delta_ms": None, "runtime_delta_pct": None}
    try:
        rb_f = float(rb)
        ra_f = float(ra)
    except (TypeError, ValueError):
        return {"runtime_delta_ms": None, "runtime_delta_pct": None}
    delta_ms = ra_f - rb_f
    if rb_f == 0:
        delta_pct = None if ra_f == 0 else float("inf")
    else:
        delta_pct = (delta_ms / rb_f) * 100.0
    return {"runtime_delta_ms": delta_ms, "runtime_delta_pct": delta_pct}


def build_comparison_report(
    *,
    before_profile: dict[str, Any],
    after_profile: dict[str, Any],
    before_mode: str = "full",
    after_mode: str = "incremental",
) -> dict[str, Any]:
    """
    Build before/after comparison report for trio only.

    Args:
        before_profile: Profile from run_baseline_profile with apsp_mode=full.
        after_profile: Profile from run_baseline_profile with apsp_mode=incremental.
        before_mode: Label for before (e.g. "full").
        after_mode: Label for after (e.g. "incremental").

    Returns:
        Report dict with metadata, per_scale comparison, failure_breakdown, interpretation.
    """
    before_per_scale = filter_trio_only(before_profile.get("per_scale", {}))
    after_per_scale = filter_trio_only(after_profile.get("per_scale", {}))

    scales = sorted(before_per_scale.keys(), key=int)
    comparison: dict[str, Any] = {}
    failure_breakdown: list[dict[str, Any]] = []

    for n_str in scales:
        b_scale = before_per_scale.get(n_str, {})
        a_scale = after_per_scale.get(n_str, {})
        b_runs = b_scale.get("runs", {})
        a_runs = a_scale.get("runs", {})

        scale_rows: dict[str, Any] = {}
        for algo_id in TRIO_ALGORITHMS:
            b_data = b_runs.get(algo_id, {})
            a_data = a_runs.get(algo_id, {})

            b_metrics = _extract_run_metrics(b_data)
            a_metrics = _extract_run_metrics(a_data)

            b_failed = bool(b_data.get("errors")) or b_metrics.get("status") is False
            a_failed = bool(a_data.get("errors")) or a_metrics.get("status") is False

            if b_failed:
                failure_breakdown.append({
                    "scale": int(n_str),
                    "algorithm_id": algo_id,
                    "mode": before_mode,
                    "error": b_data.get("errors", ["unknown"])[0] if b_data.get("errors") else "status=False",
                })
            if a_failed:
                failure_breakdown.append({
                    "scale": int(n_str),
                    "algorithm_id": algo_id,
                    "mode": after_mode,
                    "error": a_data.get("errors", ["unknown"])[0] if a_data.get("errors") else "status=False",
                })

            deltas = _compute_deltas(b_data, a_data) if not (b_failed or a_failed) else {
                "runtime_delta_ms": None,
                "runtime_delta_pct": None,
            }

            scale_rows[algo_id] = {
                "before": {
                    "runtime_ms_median": b_metrics.get("runtime_ms_median"),
                    "max_stretch": b_metrics.get("max_stretch"),
                    "status": b_metrics.get("status"),
                    "verification_policy": before_profile.get("verification_policy", "strict"),
                },
                "after": {
                    "runtime_ms_median": a_metrics.get("runtime_ms_median"),
                    "max_stretch": a_metrics.get("max_stretch"),
                    "status": a_metrics.get("status"),
                    "verification_policy": after_profile.get("verification_policy", "strict"),
                },
                "runtime_delta_ms": deltas["runtime_delta_ms"],
                "runtime_delta_pct": deltas["runtime_delta_pct"],
                "correctness_match": (
                    b_metrics.get("max_stretch") == a_metrics.get("max_stretch")
                    and b_metrics.get("status") == a_metrics.get("status")
                    if not (b_failed or a_failed) else None
                ),
            }
        comparison[n_str] = {"n": int(n_str), "runs": scale_rows}

    # Interpretation and rollout recommendation
    all_ok = len(failure_breakdown) == 0
    any_correctness_mismatch = False
    for n_str, scale_data in comparison.items():
        for algo_id, row in scale_data["runs"].items():
            if row.get("correctness_match") is False:
                any_correctness_mismatch = True
                break

    interpretation = []
    if all_ok and not any_correctness_mismatch:
        interpretation.append("All trio runs passed with strict verification.")
        interpretation.append("Correctness: before and after produce identical max_stretch and status.")
        interpretation.append("Rollout recommendation: incremental APSP mode is safe for trio at these scales.")
    elif failure_breakdown:
        interpretation.append(f"Failures: {len(failure_breakdown)} run(s) failed.")
        interpretation.append("Rollout recommendation: investigate failures before rollout.")
    elif any_correctness_mismatch:
        interpretation.append("Correctness mismatch: before and after differ in max_stretch or status.")
        interpretation.append("Rollout recommendation: do not rollout until correctness is verified.")

    report = {
        "report_type": "before_after_comparison",
        "algorithms": list(TRIO_ALGORITHMS),
        "trio_only": True,
        "before_mode": before_mode,
        "after_mode": after_mode,
        "verification_policy": "strict",
        "scales": [int(s) for s in scales],
        "seed": before_profile.get("seed"),
        "repeats": before_profile.get("repeats"),
        "t": before_profile.get("t"),
        "problem_type": before_profile.get("problem_type", "t"),
        "reproducibility": {
            "deterministic_seeds": True,
            "repeat_strategy": "seed + scale + repeat_idx per run",
        },
        "per_scale": comparison,
        "failure_breakdown": failure_breakdown,
        "interpretation": interpretation,
        "rollout_recommendation": interpretation[-1] if interpretation else "Unknown",
    }
    return report


def write_comparison_report(
    out_dir: Path,
    report: dict[str, Any],
    *,
    timestamp: datetime | None = None,
) -> Path:
    """
    Write comparison_report.json and comparison_report.md to out_dir.

    Args:
        out_dir: Directory to write artifacts.
        report: Report dict from build_comparison_report.
        timestamp: Optional timestamp for header.

    Returns:
        out_dir path.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ts = timestamp or datetime.now()
    ts_str = ts.strftime("%Y-%m-%dT%H:%M:%S")

    json_path = out_dir / "comparison_report.json"
    with open(json_path, "w", encoding="utf-8") as f:
        safe_dump(report, f, indent=2)

    md_lines = [
        "# SP-59 Before/After Comparison Report (Trio)",
        "",
        f"**Generated:** {ts_str}",
        f"**Algorithms (trio only):** {report.get('algorithms', [])}",
        f"**Before mode:** {report.get('before_mode')}",
        f"**After mode:** {report.get('after_mode')}",
        f"**Verification policy:** {report.get('verification_policy')}",
        f"**Scales:** {report.get('scales', [])}",
        f"**Seed:** {report.get('seed')}",
        f"**t:** {report.get('t')}",
        "",
        "## Per-scale comparison",
        "",
    ]

    for n_str, scale_data in report.get("per_scale", {}).items():
        md_lines.append(f"### n = {n_str}")
        md_lines.append("")
        for algo_id, row in scale_data.get("runs", {}).items():
            md_lines.append(f"- **{algo_id}**")
            b = row.get("before", {})
            a = row.get("after", {})
            md_lines.append(f"  - Before: runtime_ms={b.get('runtime_ms_median')}, max_stretch={b.get('max_stretch')}, status={b.get('status')}")
            md_lines.append(f"  - After:  runtime_ms={a.get('runtime_ms_median')}, max_stretch={a.get('max_stretch')}, status={a.get('status')}")
            md_lines.append(f"  - runtime_delta_ms={row.get('runtime_delta_ms')}, runtime_delta_pct={row.get('runtime_delta_pct')}")
            md_lines.append(f"  - correctness_match={row.get('correctness_match')}")
            md_lines.append("")
        md_lines.append("")

    if report.get("failure_breakdown"):
        md_lines.append("## Failure breakdown")
        md_lines.append("")
        for fb in report["failure_breakdown"]:
            md_lines.append(f"- n={fb.get('scale')} {fb.get('algorithm_id')} ({fb.get('mode')}): {fb.get('error')}")
        md_lines.append("")

    md_lines.append("## Interpretation")
    md_lines.append("")
    for line in report.get("interpretation", []):
        md_lines.append(f"- {line}")
    md_lines.append("")
    md_lines.append(f"**Rollout recommendation:** {report.get('rollout_recommendation', 'N/A')}")
    md_lines.append("")

    md_path = out_dir / "comparison_report.md"
    with open(md_path, "w", encoding="utf-8") as f:
        f.write("\n".join(md_lines))

    return out_dir
