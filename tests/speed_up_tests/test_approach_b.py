"""
Approach B tests: CuBlock polynomial fit limited to input sample only (j=1:1).

Since Approach B produces mathematically identical output to the original
(only the wasted P-column fits are skipped), it validates against
expected_output.csv with the same atol=1e-4 as test_single_worker_matches_reference.
"""

import pytest

from .conftest import (
    compare_to_expected,
    load_fixtures,
    measure_wall_time,
    octave_required,
    run_pipeline,
)


@octave_required
@pytest.mark.timeout(300)
def test_b_output_identical_to_baseline():
    """
    Approach B output must match expected_output.csv within atol=1e-4.

    Since expected_output was generated with the original CuBlock (all-columns loop),
    matching within 1e-4 confirms B produces mathematically identical results.
    """
    input_df, p_df, q_df, expected_df = load_fixtures()
    result = run_pipeline(input_df, p_df, q_df, n_workers=1, random_seed=42)
    metrics = compare_to_expected(result, expected_df)

    print(f"\nApproach B faithfulness metrics:")
    print(f"  max_abs_diff  = {metrics['max_abs_diff']:.8f}")
    print(f"  mean_abs_diff = {metrics['mean_abs_diff']:.8f}")
    print(f"  mean_rel_diff = {metrics['mean_rel_diff']:.6%}")
    print(f"  %within_1pct  = {metrics['pct_within_1pct']:.1f}%")

    assert metrics["max_abs_diff"] < 1e-4, (
        f"Approach B max_abs_diff={metrics['max_abs_diff']:.2e} exceeds atol=1e-4. "
        f"The j=1:1 fix must not change output for the input sample column."
    )


@octave_required
@pytest.mark.timeout(300)
def test_b_timing():
    """Log wall-clock time for Approach B (post-fix). Informational — no assertion."""
    input_df, p_df, q_df, _ = load_fixtures()

    with measure_wall_time() as t:
        run_pipeline(input_df, p_df, q_df, n_workers=1, random_seed=42)

    elapsed = t[0]
    n_samples = input_df.shape[0]
    per_sample = elapsed / n_samples
    print(f"\nApproach B timing (10 samples, 1 worker, P0_small NP=39):")
    print(f"  Total wall time : {elapsed:.1f} s")
    print(f"  Per-sample time : {per_sample:.2f} s/sample")
    print(f"  (Compare against NBGPL570 baseline ~155 s/sample to estimate full speedup)")
