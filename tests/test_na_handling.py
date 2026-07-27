"""
Unit tests for shambhala/na_handling.py.
"""

import pathlib
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from shambhala.na_handling import (
    NAMask,
    apply_na_strategy,
    detect_na_genes,
    restore_na_genes,
)


def _make_df(n_samples: int = 6, n_genes: int = 10, seed: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    data = rng.random((n_samples, n_genes)) * 100
    return pd.DataFrame(
        data,
        index=[f"S{i}" for i in range(n_samples)],
        columns=[f"G{j}" for j in range(n_genes)],
    )


def _inject_nans(df: pd.DataFrame, genes: list[str], sample_idx: int = 0) -> pd.DataFrame:
    df = df.copy()
    for g in genes:
        df.loc[df.index[sample_idx], g] = np.nan
    return df


def test_detect_na_genes_correct():
    df = _make_df()
    na_genes = ["G2", "G5", "G8"]
    df = _inject_nans(df, na_genes)
    detected = detect_na_genes(df)
    assert set(detected) == set(na_genes)


def test_drop_removes_na_genes():
    df = _make_df()
    df = _inject_nans(df, ["G1", "G3"])
    clean_df, mask = apply_na_strategy(df, strategy="drop")
    assert clean_df.isna().sum().sum() == 0
    assert "G1" not in clean_df.columns
    assert "G3" not in clean_df.columns
    assert mask.strategy == "drop"
    assert set(mask.dropped_genes) == {"G1", "G3"}


def test_knn_fills_nas():
    df = _make_df()
    df = _inject_nans(df, ["G2", "G7"])
    clean_df, mask = apply_na_strategy(df, strategy="knn", knn_k=3)
    assert clean_df.isna().sum().sum() == 0
    assert mask.strategy == "knn"
    assert "G2" in clean_df.columns
    assert "G7" in clean_df.columns


def test_restore_preserves_original_order():
    df = _make_df()
    original_cols = list(df.columns)
    df = _inject_nans(df, ["G4"])
    clean_df, mask = apply_na_strategy(df, strategy="drop")
    restored = restore_na_genes(clean_df, mask)
    assert list(restored.columns) == original_cols


def test_restore_inserts_nan_correct_columns():
    df = _make_df()
    df = _inject_nans(df, ["G0", "G9"])
    clean_df, mask = apply_na_strategy(df, strategy="drop")
    restored = restore_na_genes(clean_df, mask)
    assert restored["G0"].isna().all()
    assert restored["G9"].isna().all()
    # Non-dropped genes should have no NaN
    for col in restored.columns:
        if col not in {"G0", "G9"}:
            assert not restored[col].isna().any(), f"Unexpected NaN in column {col}"


def test_all_na_gene_handled_gracefully():
    df = _make_df()
    # G5 has 100% NaN — exceeds max_na_frac=0.20, so it should be dropped
    df["G5"] = np.nan
    clean_df, mask = apply_na_strategy(df, strategy="knn", knn_k=3, max_na_frac=0.20)
    assert "G5" not in clean_df.columns
    assert "G5" in mask.dropped_genes
    assert clean_df.isna().sum().sum() == 0
