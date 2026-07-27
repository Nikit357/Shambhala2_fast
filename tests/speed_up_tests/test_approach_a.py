"""
Approach A tests: precompute QN reference from P; apply Python-side QN per sample.

--precompute-qn-reference replaces Octave's quantilenorm on G×(1+NP) with a
Python C-extension call on G×1 (each sample against precomputed P rank means).
--synthetic-cublock-p additionally passes 1 synthetic P column to Octave CuBlock.

Approach A+synthetic is documented as exploratory only (~15% mean_rel_diff expected).
"""

import pytest

from .conftest import (
    assert_mean_rel_diff,
    compare_to_expected,
    load_fixtures,
    measure_wall_time,
    octave_required,
    run_pipeline,
)

try:
    import qnorm  # noqa: F401
    _QNORM_AVAILABLE = True
except ImportError:
    _QNORM_AVAILABLE = False

qnorm_required = pytest.mark.skipif(
    not _QNORM_AVAILABLE,
    reason="qnorm library not installed (pip install qnorm)",
)


@octave_required
@qnorm_required
@pytest.mark.timeout(3600)
def test_a_qn_only():
    """Approach A (--precompute-qn-reference only): mean_rel_diff < 5%."""
    input_df, p_df, q_df, expected_df = load_fixtures()
    result = run_pipeline(
        input_df, p_df, q_df,
        n_workers=1,
        random_seed=42,
        precompute_qn_reference=True,
        synthetic_cublock_p=False,
    )
    metrics = compare_to_expected(result, expected_df)
    assert_mean_rel_diff(metrics, warn_threshold=0.01, fail_threshold=0.05, label="Approach A (QN only)")


@octave_required
@qnorm_required
@pytest.mark.timeout(3600)
def test_a_synthetic_cublock_p():
    """
    Approach A + --synthetic-cublock-p: mean_rel_diff < 20%.

    This approach is EXPLORATORY ONLY. ~15% mean_rel_diff is expected because
    replacing 39 P columns with a single synthetic column collapses CuBlock's
    gene k-means from 40D to 2D, introducing large systematic deviation.
    Do not use this flag for production benchmarking.
    """
    input_df, p_df, q_df, expected_df = load_fixtures()
    result = run_pipeline(
        input_df, p_df, q_df,
        n_workers=1,
        random_seed=42,
        precompute_qn_reference=True,
        synthetic_cublock_p=True,
    )
    metrics = compare_to_expected(result, expected_df)
    assert_mean_rel_diff(
        metrics,
        warn_threshold=0.10,
        fail_threshold=0.20,
        label="Approach A + synthetic P (EXPLORATORY)",
    )


@octave_required
@qnorm_required
@pytest.mark.timeout(3600)
@pytest.mark.slow
def test_a_timing_qn_only():
    """Log wall-clock time for Approach A (QN only). Informational — no assertion."""
    input_df, p_df, q_df, _ = load_fixtures()

    with measure_wall_time() as t_baseline:
        run_pipeline(input_df, p_df, q_df, n_workers=1, random_seed=42)
    with measure_wall_time() as t_a:
        run_pipeline(
            input_df, p_df, q_df, n_workers=1, random_seed=42,
            precompute_qn_reference=True,
        )

    speedup = t_baseline[0] / t_a[0] if t_a[0] > 0 else float("inf")
    print(f"\nApproach A (QN only) timing:")
    print(f"  Baseline: {t_baseline[0]:.1f} s  |  Approach A: {t_a[0]:.1f} s  |  Speedup: {speedup:.2f}×")


@octave_required
@qnorm_required
@pytest.mark.timeout(3600)
@pytest.mark.slow
def test_a_timing_combined():
    """Log wall-clock time for Approach A + synthetic CuBlock P. Informational."""
    input_df, p_df, q_df, _ = load_fixtures()

    with measure_wall_time() as t_baseline:
        run_pipeline(input_df, p_df, q_df, n_workers=1, random_seed=42)
    with measure_wall_time() as t_a:
        run_pipeline(
            input_df, p_df, q_df, n_workers=1, random_seed=42,
            precompute_qn_reference=True,
            synthetic_cublock_p=True,
        )

    speedup = t_baseline[0] / t_a[0] if t_a[0] > 0 else float("inf")
    print(f"\nApproach A + synthetic_cublock_p timing:")
    print(f"  Baseline: {t_baseline[0]:.1f} s  |  Combined: {t_a[0]:.1f} s  |  Speedup: {speedup:.2f}×")
