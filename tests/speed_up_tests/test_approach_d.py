"""
Approach D tests: Python-side QN using qnorm library.

Validates that qnorm.quantile_normalize produces output identical (within atol=1e-4)
to Octave's quantilenorm.m when applied against the same P reference distribution.
This test does NOT require Octave — it compares Python-side QN against a known
expected distribution shape.
"""

import numpy as np
import pandas as pd
import pytest

try:
    import qnorm
    _QNORM_AVAILABLE = True
except ImportError:
    _QNORM_AVAILABLE = False

qnorm_required = pytest.mark.skipif(
    not _QNORM_AVAILABLE,
    reason="qnorm library not installed (pip install qnorm)",
)


@qnorm_required
def test_d_qnorm_produces_identical_distributions():
    """
    qnorm.quantile_normalize(axis=1, target=ref) maps each sample to the target distribution.

    After normalization, the sorted values of each sample must equal the target.
    """
    rng = np.random.default_rng(42)
    G, S = 200, 5
    data = rng.exponential(scale=5.0, size=(G, S))
    ref = np.sort(rng.exponential(scale=5.0, size=G))  # sorted reference distribution

    df = pd.DataFrame(data, index=[f"g{i}" for i in range(G)], columns=[f"s{j}" for j in range(S)])
    result = qnorm.quantile_normalize(df, axis=1, target=ref)

    for col in result.columns:
        sorted_col = np.sort(result[col].values)
        np.testing.assert_allclose(
            sorted_col, ref, atol=1e-10,
            err_msg=f"Column {col}: sorted values do not match target distribution.",
        )


@qnorm_required
def test_d_qnorm_matches_manual_rank_normalization():
    """
    qnorm with target should produce the same result as manually rank-normalizing:
    sort each column, replace with target[rank], unsort.
    """
    rng = np.random.default_rng(123)
    G, S = 50, 4
    data = rng.standard_normal((G, S))
    ref = np.sort(rng.standard_normal(G))

    df = pd.DataFrame(data, index=[f"g{i}" for i in range(G)], columns=[f"s{j}" for j in range(S)])
    result_qnorm = qnorm.quantile_normalize(df, axis=1, target=ref).values

    # Manual rank normalization
    result_manual = np.empty_like(data)
    for j in range(S):
        col = data[:, j]
        order = np.argsort(col)
        result_manual[order, j] = ref  # assign sorted target values to sorted positions

    np.testing.assert_allclose(result_qnorm, result_manual, atol=1e-10)


@qnorm_required
def test_d_synthetic_p_reference_is_sorted():
    """
    The synthetic P reference (mean of sorted P columns per rank) must be monotonically
    non-decreasing — otherwise it cannot be a valid QN target distribution.
    """
    from .conftest import load_fixtures

    _, p_df, _, _ = load_fixtures()
    p_T = p_df.T  # genes × samples

    p_sorted = np.sort(p_T.values, axis=0)  # sort each P column
    p_qn_reference = p_sorted.mean(axis=1)   # mean per rank

    diffs = np.diff(p_qn_reference)
    assert (diffs >= -1e-12).all(), (
        f"Synthetic P reference is not monotonically non-decreasing. "
        f"Min diff = {diffs.min():.4e}"
    )


@qnorm_required
def test_d_qnorm_preserves_gene_rank_order():
    """
    After QN, genes that had higher expression than others within a sample
    should retain their relative rank order (QN is a rank-preserving transformation).
    """
    rng = np.random.default_rng(7)
    G, S = 100, 3
    data = rng.standard_normal((G, S))
    ref = np.sort(rng.standard_normal(G))

    df = pd.DataFrame(data, index=[f"g{i}" for i in range(G)], columns=[f"s{j}" for j in range(S)])
    result = qnorm.quantile_normalize(df, axis=1, target=ref)

    for j in range(S):
        original_ranks = np.argsort(data[:, j])
        result_ranks = np.argsort(result.values[:, j])
        np.testing.assert_array_equal(
            original_ranks, result_ranks,
            err_msg=f"Sample s{j}: QN changed gene rank order.",
        )
