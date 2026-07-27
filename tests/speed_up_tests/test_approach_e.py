"""
Tests for Approach E: precompute k-means gene clusters from P once; pass to CuBlock_fixed.

Test matrix:
  E-1  test_e_clusters_stable         — same seed → identical cluster labels (pure Python)
  E-2  test_e_clusters_shape          — labels are 1-indexed int32, length == n_genes (pure Python)
  E-3  test_e_output_within_tolerance — mean_rel_diff < 5% vs expected_output.csv (Octave)
  E-4  test_e_timing                  — log wall time; no hard assertion
"""

import numpy as np
import pytest

from tests.speed_up_tests.conftest import (
    assert_mean_rel_diff,
    compare_to_expected,
    load_fixtures,
    measure_wall_time,
    octave_required,
    run_pipeline,
)
from shambhala.octave_bridge import precompute_clusters_from_p


def test_e_clusters_stable():
    _, p_df, _, _ = load_fixtures()
    p_T = p_df.T
    labels_a = precompute_clusters_from_p(p_T=p_T, k=5, random_seed=7)
    labels_b = precompute_clusters_from_p(p_T=p_T, k=5, random_seed=7)
    np.testing.assert_array_equal(labels_a, labels_b)


def test_e_clusters_shape():
    _, p_df, _, _ = load_fixtures()
    p_T = p_df.T
    n_genes = p_T.shape[0]
    labels = precompute_clusters_from_p(p_T=p_T, k=5, random_seed=42)
    assert labels.shape == (n_genes,)
    assert labels.dtype == np.int32
    assert labels.min() >= 1
    assert labels.max() <= 5


@octave_required
@pytest.mark.timeout(3600)
def test_e_output_within_tolerance():
    input_df, p_df, q_df, expected_df = load_fixtures()
    result = run_pipeline(
        input_df=input_df,
        p_df=p_df,
        q_df=q_df,
        n_workers=1,
        precompute_cublock_clusters=True,
        random_seed=42,
    )
    metrics = compare_to_expected(result, expected_df)
    assert_mean_rel_diff(metrics, warn_threshold=0.01, fail_threshold=0.05, label="Approach E")


@octave_required
@pytest.mark.timeout(3600)
@pytest.mark.slow
def test_e_timing():
    input_df, p_df, q_df, _ = load_fixtures()
    with measure_wall_time() as t_base:
        run_pipeline(input_df=input_df, p_df=p_df, q_df=q_df, n_workers=1, random_seed=42)
    with measure_wall_time() as t_e:
        run_pipeline(
            input_df=input_df,
            p_df=p_df,
            q_df=q_df,
            n_workers=1,
            precompute_cublock_clusters=True,
            random_seed=42,
        )
    print(
        f"\nApproach E timing: baseline={t_base[0]:.1f}s  with_clusters={t_e[0]:.1f}s"
        f"  ratio={t_base[0] / max(t_e[0], 0.001):.2f}x"
    )
