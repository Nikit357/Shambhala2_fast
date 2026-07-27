"""
Unit tests for shambhala/q_rescale.py.
"""

import pathlib
import sys
import warnings

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from shambhala.q_rescale import compute_q_statistics, rescale


def _make_q_df(n_samples: int = 5, n_genes: int = 8, seed: int = 1) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    data = rng.random((n_samples, n_genes)) * 1000 + 1
    return pd.DataFrame(
        data,
        index=[f"QS{i}" for i in range(n_samples)],
        columns=[f"G{j}" for j in range(n_genes)],
    )


def test_rescale_known_values():
    # Hand-computed expected output for a 2-gene, 1-sample toy case.
    # Q gene A: values [2.0, 4.0]  → log([2+1e-6, 4+1e-6])
    # rm_A ≈ mean(log([2, 4])) = (ln(2)+ln(4))/2 ≈ 1.03972
    # rs_A ≈ std(log([2, 4])) with ddof=1 = ln(2) ≈ 0.69315
    q_df = pd.DataFrame({"A": [2.0, 4.0]}, index=["QS0", "QS1"])
    rm, rs = compute_q_statistics(q_df, q_pseudocount=0.0)

    # octave_output gene A = 3.0
    octave_out = pd.DataFrame({"S0": [3.0]}, index=["A"])
    result = rescale(octave_out, rm, rs)

    expected_log = rm["A"] + rs["A"] * np.log(3.0 + 1)
    expected = np.exp(expected_log)
    assert abs(result.loc["A", "S0"] - expected) < 1e-10


def test_compute_q_statistics_shape():
    q_df = _make_q_df()
    rm, rs = compute_q_statistics(q_df)
    assert isinstance(rm, pd.Series)
    assert isinstance(rs, pd.Series)
    assert list(rm.index) == list(q_df.columns)
    assert list(rs.index) == list(q_df.columns)


def test_pseudocount_prevents_log_zero():
    q_df = pd.DataFrame({"A": [0.0, 5.0], "B": [3.0, 4.0]}, index=["QS0", "QS1"])
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        rm, rs = compute_q_statistics(q_df, q_pseudocount=1e-6)
        assert len(w) == 0, "Unexpected warning raised with default pseudocount"
    assert "A" in rm.index
    assert not np.isnan(rm["A"])


def test_pseudocount_zero_excludes_with_warning():
    q_df = pd.DataFrame({"A": [0.0, 5.0], "B": [3.0, 4.0]}, index=["QS0", "QS1"])
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        rm, rs = compute_q_statistics(q_df, q_pseudocount=0.0)
        assert len(w) == 1
        assert "A" in str(w[0].message)
    assert "A" not in rm.index
    assert "B" in rm.index


def test_rescale_raises_on_negative_q():
    q_df = pd.DataFrame({"A": [-1.0, 5.0]}, index=["QS0", "QS1"])
    with pytest.raises(ValueError, match="negative values"):
        compute_q_statistics(q_df)


def test_rescale_gene_intersection():
    q_df = _make_q_df(n_genes=4)
    rm, rs = compute_q_statistics(q_df)
    # octave_output has 5 genes, but only 4 are in Q
    octave_out = pd.DataFrame(
        np.ones((5, 2)),
        index=["G0", "G1", "G2", "G3", "EXTRA"],
        columns=["S0", "S1"],
    )
    result = rescale(octave_out, rm, rs)
    assert "EXTRA" not in result.index
    assert result.shape[0] == 4
