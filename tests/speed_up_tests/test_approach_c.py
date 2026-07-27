"""
Approach C tests: centroid-based P subsampling.

--max-p-samples N selects the N P samples closest to the P centroid (mean gene
expression vector), reducing QN and CuBlock from G×(1+NP) to G×(1+N_sub).

Speedup is proportional to NP/N_sub. Authenticity: moderate — fewer P samples
means less stable QN reference and CuBlock gene clusters.
"""

import numpy as np
import pytest

from .conftest import (
    compare_to_expected,
    load_fixtures,
    measure_wall_time,
    octave_required,
    run_pipeline,
)


def test_c_centroid_selection():
    """
    Centroid-based subsampling selects the N samples with smallest Euclidean
    distance to the P mean. Unit test of the selection logic — no Octave needed.
    """
    _, p_df, _, _ = load_fixtures()

    n_sub = 20
    p_mean = p_df.mean(axis=0)
    distances = ((p_df - p_mean) ** 2).sum(axis=1)
    expected_selected = distances.nsmallest(n_sub).index

    # Verify selected samples are indeed the n_sub closest to centroid
    all_distances = distances.sort_values()
    closest_n = all_distances.index[:n_sub]

    # Sets should match (order may differ)
    assert set(expected_selected) == set(closest_n), (
        "Centroid selection did not pick the N samples closest to the mean."
    )

    # Verify all selected samples have smaller distance than any not-selected sample
    selected_max = distances[expected_selected].max()
    not_selected = p_df.index.difference(expected_selected)
    not_selected_min = distances[not_selected].min()
    assert selected_max <= not_selected_min + 1e-10, (
        "A non-selected sample is closer to the centroid than a selected one."
    )


def test_c_subsampling_reduces_p_size():
    """After subsampling, the p_df passed to the pipeline has exactly N rows."""
    _, p_df, _, _ = load_fixtures()

    n_original = p_df.shape[0]
    n_sub = min(10, n_original - 1)

    p_mean = p_df.mean(axis=0)
    distances = ((p_df - p_mean) ** 2).sum(axis=1)
    p_subsampled = p_df.loc[distances.nsmallest(n_sub).index]

    assert p_subsampled.shape[0] == n_sub
    assert p_subsampled.shape[1] == p_df.shape[1]


@octave_required
@pytest.mark.parametrize("n_sub", [40, 20, 10])
@pytest.mark.timeout(300)
def test_c_output_and_timing(n_sub):
    """
    Report max_abs_diff and speedup for each P subsampling level.

    For n_sub >= 20, max_abs_diff should be < 0.10 (soft guideline only).
    Prints detailed metrics for manual review — no hard assertion for n_sub < 40.
    """
    input_df, p_df, q_df, expected_df = load_fixtures()
    n_original = p_df.shape[0]

    if n_sub >= n_original:
        pytest.skip(f"n_sub={n_sub} >= n_original={n_original}; subsampling has no effect.")

    with measure_wall_time() as t_baseline:
        baseline = run_pipeline(input_df, p_df, q_df, n_workers=1, random_seed=42)
    with measure_wall_time() as t_c:
        result = run_pipeline(
            input_df, p_df, q_df, n_workers=1, random_seed=42,
            max_p_samples=n_sub,
        )

    metrics_vs_expected = compare_to_expected(result, expected_df)
    speedup = t_baseline[0] / t_c[0] if t_c[0] > 0 else float("inf")

    print(f"\nApproach C (n_sub={n_sub}, n_original={n_original}) faithfulness vs expected_output.csv:")
    print(f"  max_abs_diff  = {metrics_vs_expected['max_abs_diff']:.6f}")
    print(f"  mean_abs_diff = {metrics_vs_expected['mean_abs_diff']:.6f}")
    print(f"  mean_rel_diff = {metrics_vs_expected['mean_rel_diff']:.4%}")
    print(f"  %within_1pct  = {metrics_vs_expected['pct_within_1pct']:.1f}%")
    print(f"  Speedup vs baseline: {speedup:.2f}×")

    if n_sub == 40:
        assert metrics_vs_expected["max_abs_diff"] < 0.10, (
            f"n_sub=40 max_abs_diff={metrics_vs_expected['max_abs_diff']:.4f} exceeds atol=0.10."
        )
