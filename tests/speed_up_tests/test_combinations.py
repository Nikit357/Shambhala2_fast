"""
Phase 6 — Combination tests + summary table.

Runs the 10-sample fixture through selected combinations of speed-up flags.
All tests require Octave; pure-Python checks are deferred to per-approach files.

The session-level baseline timing is cached in `conftest.py` via the
`baseline_result` fixture so all combination tests compare against the same run.

Combination IDs match the summary table in Section 10 of the plan document.
"""

import pytest

from tests.speed_up_tests.conftest import (
    assert_mean_rel_diff,
    compare_to_expected,
    load_fixtures,
    measure_wall_time,
    octave_required,
    print_expression_comparison,
    run_pipeline,
)

_COMBO_RESULTS: dict[str, dict] = {}

_COMBOS: list = [
    pytest.param(
        "baseline", {}, (0.01, 0.05), "Baseline: Approach B only (always on)", id="baseline"
    ),
    pytest.param(
        "A_only", {"precompute_qn_reference": True}, (0.01, 0.05), "A: Python QN reference",
        id="A_only"
    ),
    pytest.param(
        "A_D",
        {"precompute_qn_reference": True, "synthetic_cublock_p": True},
        (0.10, 0.20),
        "A+D: Python QN + synthetic 1-col P for CuBlock (EXPLORATORY)",
        id="A_D",
    ),
    pytest.param(
        "C_40", {"max_p_samples": 40}, (0.01, 0.05), "C: centroid P subsample to 40", id="C_40"
    ),
    pytest.param(
        "C_20", {"max_p_samples": 20}, (0.01, 0.05), "C: centroid P subsample to 20", id="C_20"
    ),
    pytest.param(
        "E",
        {"precompute_cublock_clusters": True},
        (0.01, 0.05),
        "E: precomputed k-means clusters",
        id="E",
    ),
    pytest.param(
        "A_E",
        {"precompute_qn_reference": True, "precompute_cublock_clusters": True},
        (0.01, 0.05),
        "A+E: Python QN + precomputed clusters",
        id="A_E",
    ),
    pytest.param(
        "G",
        {"python_cublock": True},
        (None, None),
        "G: Python CuBlock",
        marks=pytest.mark.xfail(
            strict=False,
            reason="python_cublock=True triggers ProcessPoolExecutor deadlock in harmonize_parallel",
        ),
        id="G",
    ),
    pytest.param(
        "A_D_E",
        {
            "precompute_qn_reference": True,
            "synthetic_cublock_p": True,
            "precompute_cublock_clusters": True,
        },
        (0.10, 0.20),
        "A+D+E: Python QN + synthetic P + precomputed clusters (EXPLORATORY)",
        id="A_D_E",
    ),
    pytest.param(
        "A_D_G",
        {"precompute_qn_reference": True, "python_cublock": True, "synthetic_cublock_p": True},
        (None, None),
        "A+D+G: Python QN + synthetic P + Python CuBlock",
        marks=pytest.mark.xfail(
            strict=False,
            reason="python_cublock=True triggers ProcessPoolExecutor deadlock in harmonize_parallel",
        ),
        id="A_D_G",
    ),
]


@pytest.fixture(scope="session")
def baseline_result():
    input_df, p_df, q_df, expected_df = load_fixtures()
    with measure_wall_time() as t:
        result = run_pipeline(
            input_df=input_df, p_df=p_df, q_df=q_df, n_workers=1, random_seed=42
        )
    return result, expected_df, t[0]


@octave_required
@pytest.mark.parametrize("combo_id,kwargs,thresholds,desc", _COMBOS)
def test_combination(combo_id, kwargs, thresholds, desc, baseline_result):
    baseline_df, expected_df, baseline_time = baseline_result
    input_df, p_df, q_df, _ = load_fixtures()
    with measure_wall_time() as t:
        result = run_pipeline(
            input_df=input_df, p_df=p_df, q_df=q_df, n_workers=1, random_seed=42, **kwargs
        )
    metrics = compare_to_expected(result, expected_df)
    speedup = baseline_time / max(t[0], 1e-3)
    warn_threshold, fail_threshold = thresholds

    mrd = metrics["mean_rel_diff"]
    if mrd < warn_threshold:
        status = "PASS"
    elif mrd < fail_threshold:
        status = "WARN"
    else:
        status = "FAIL"

    _COMBO_RESULTS[combo_id] = {
        "desc": desc,
        "mean_rel_diff": mrd,
        "max_abs_diff": metrics["max_abs_diff"],
        "mean_abs_diff": metrics["mean_abs_diff"],
        "pct_within_1pct": metrics["pct_within_1pct"],
        "time_s": t[0],
        "speedup": speedup,
        "status": status,
        "warn": warn_threshold,
        "fail": fail_threshold,
    }

    print(
        f"\n[{combo_id}] {desc}\n"
        f"  mean_rel_diff={mrd:.4%}  max_abs_diff={metrics['max_abs_diff']:.4f}"
        f"  mean_abs_diff={metrics['mean_abs_diff']:.4f}"
        f"  pct_within_1pct={metrics['pct_within_1pct']:.1f}%"
        f"  time={t[0]:.1f}s  speedup={speedup:.2f}x  status={status}"
    )

    print_expression_comparison(result, expected_df, label=combo_id, n_genes=8, n_samples=3)

    assert_mean_rel_diff(
        metrics,
        warn_threshold=warn_threshold,
        fail_threshold=fail_threshold,
        label=combo_id,
    )


def _combo_id(entry: object) -> str:
    return str(getattr(entry, "id", ""))  # all entries are pytest.param


@octave_required
def test_z_print_summary_table():
    divider = "─" * 97
    print(f"\n{divider}")
    print(
        f"{'COMBO':<12}  {'mean_rel_diff':>14}  {'max_abs_diff':>13}  "
        f"{'mean_abs_diff':>13}  {'%within_1pct':>12}  {'time(s)':>7}  {'speedup':>8}  status"
    )
    print(divider)

    for entry in _COMBOS:
        cid = _combo_id(entry)
        if cid not in _COMBO_RESULTS:
            print(
                f"{cid:<12}  {'—':>14}  {'—':>13}  {'—':>13}  {'—':>12}  {'—':>7}  {'—':>8}  XFAIL / TIMEOUT"
            )
            continue
        r = _COMBO_RESULTS[cid]
        print(
            f"{cid:<12}  {r['mean_rel_diff']:>13.4%}  {r['max_abs_diff']:>13.4f}  "
            f"{r['mean_abs_diff']:>13.4f}  {r['pct_within_1pct']:>11.1f}%  "
            f"{r['time_s']:>7.1f}  {r['speedup']:>7.2f}x  {r['status']}"
        )

    print(divider)
    print("  Legend: PASS < 1% | WARN 1%–5% | FAIL > 5%  (exploratory combos: WARN < 10% | FAIL > 20%)")
    print()
