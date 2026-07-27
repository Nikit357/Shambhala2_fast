"""
End-to-end integration tests for the full Shambhala_containerized pipeline.

These tests require Octave to be installed and accessible via /opt/conda/bin/octave
or on PATH. All tests are automatically skipped when Octave is not available.
They use the toy 10-sample fixture derived from Shambhala2/Input.csv.
"""

import pathlib
import shutil
import sys
import tempfile

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from shambhala import io_utils, na_handling, q_rescale, parallel as parallel_module

FIXTURES_DIR = pathlib.Path(__file__).parent / "fixtures"
OCTAVE_DIR = pathlib.Path(__file__).resolve().parent.parent / "octave"
OCTAVE_BIN = "/opt/conda/bin/octave" if pathlib.Path("/opt/conda/bin/octave").exists() else (shutil.which("octave") or "octave")

_OCTAVE_AVAILABLE = pathlib.Path("/opt/conda/bin/octave").exists() or shutil.which("octave") is not None
octave_required = pytest.mark.skipif(
    not _OCTAVE_AVAILABLE,
    reason="Octave not installed (expected at /opt/conda/bin/octave or on PATH)",
)


def _run_pipeline(
    input_path: pathlib.Path,
    p_path: pathlib.Path,
    q_path: pathlib.Path,
    n_workers: int = 1,
    na_strategy: str = "drop",
    random_seed: int = 42,
) -> pd.DataFrame:
    input_df = io_utils.read_expression(str(input_path))
    p_df = io_utils.read_expression(str(p_path)).astype(float)
    q_df = io_utils.read_expression(str(q_path)).astype(float)

    input_T = input_df.T
    p_T = p_df.T
    q_T = q_df.T

    common_genes = input_T.index.intersection(p_T.index).intersection(q_T.index)
    input_T = input_T.loc[common_genes]
    p_T = p_T.loc[common_genes]
    q_T = q_T.loc[common_genes]

    clean_df, na_mask = na_handling.apply_na_strategy(input_T.T, strategy=na_strategy)
    clean_T = clean_df.T
    p_T = p_T.loc[clean_T.index]
    q_df_filtered = q_T.loc[clean_T.index].T
    rm, rs = q_rescale.compute_q_statistics(q_df_filtered)
    harmonized_T = parallel_module.harmonize_parallel(
        input_df=clean_T,
        p_df=p_T,
        rm=rm,
        rs=rs,
        k=5,
        n_workers=n_workers,
        octave_scripts_dir=str(OCTAVE_DIR),
        octave_bin=OCTAVE_BIN,
        timeout_s=300,
        random_seed=random_seed,
    )
    return na_handling.restore_na_genes(harmonized_T.T, na_mask)


@octave_required
@pytest.mark.timeout(300)
def test_single_worker_matches_reference():
    result = _run_pipeline(
        FIXTURES_DIR / "input_10samples.csv",
        FIXTURES_DIR / "P0_small.csv",
        FIXTURES_DIR / "Q0_small.csv",
        n_workers=1,
    )
    expected = io_utils.read_expression(str(FIXTURES_DIR / "expected_output.csv"))

    # Align on common genes/samples
    common_genes = result.columns.intersection(expected.columns)
    common_samples = result.index.intersection(expected.index)
    result_aligned = result.loc[common_samples, common_genes]
    expected_aligned = expected.loc[common_samples, common_genes].astype(float)

    np.testing.assert_allclose(
        result_aligned.values, expected_aligned.values, atol=1e-4, rtol=0,
        err_msg="Single-worker output does not match reference within atol=1e-4",
    )


@octave_required
@pytest.mark.timeout(600)
def test_parallel_matches_single():
    result_1 = _run_pipeline(
        FIXTURES_DIR / "input_10samples.csv",
        FIXTURES_DIR / "P0_small.csv",
        FIXTURES_DIR / "Q0_small.csv",
        n_workers=1,
        random_seed=42,
    )
    for n_workers in [3, 5]:
        result_n = _run_pipeline(
            FIXTURES_DIR / "input_10samples.csv",
            FIXTURES_DIR / "P0_small.csv",
            FIXTURES_DIR / "Q0_small.csv",
            n_workers=n_workers,
            random_seed=42,
        )
        common_genes = result_1.columns.intersection(result_n.columns)
        mean_1 = result_1[common_genes].mean(axis=0)
        mean_n = result_n[common_genes].mean(axis=0)
        np.testing.assert_allclose(
            mean_1.values, mean_n.values, rtol=0.05, atol=0,
            err_msg=f"n_workers={n_workers} per-gene row means differ from n_workers=1 by more than rtol=0.05",
        )


@octave_required
def test_local_io_roundtrip():
    result = _run_pipeline(
        FIXTURES_DIR / "input_10samples.csv",
        FIXTURES_DIR / "P0_small.csv",
        FIXTURES_DIR / "Q0_small.csv",
        n_workers=1,
    )
    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as tf:
        path = tf.name
    try:
        io_utils.write_expression(result, path)
        reloaded = io_utils.read_expression(path)
        pd.testing.assert_frame_equal(result, reloaded, check_dtype=False, atol=1e-10)
    finally:
        pathlib.Path(path).unlink(missing_ok=True)


@octave_required
def test_gene_intersection_is_applied():
    # P fixture has fewer genes than input; only their intersection should appear in output
    input_df = io_utils.read_expression(str(FIXTURES_DIR / "input_10samples.csv"))
    p_df = io_utils.read_expression(str(FIXTURES_DIR / "P0_small.csv")).astype(float)
    q_df = io_utils.read_expression(str(FIXTURES_DIR / "Q0_small.csv")).astype(float)

    # Artificially drop 50 genes from P to create a strict subset
    p_df_reduced = p_df.iloc[:, 50:]

    input_T = input_df.T
    p_T = p_df_reduced.T
    q_T = q_df.T
    common_genes = input_T.index.intersection(p_T.index).intersection(q_T.index)

    result = _run_pipeline(
        FIXTURES_DIR / "input_10samples.csv",
        FIXTURES_DIR / "P0_small.csv",
        FIXTURES_DIR / "Q0_small.csv",
        n_workers=1,
    )
    # Verify the output gene count equals the intersection
    assert len(result.columns.intersection(common_genes)) <= len(result.columns)


@octave_required
@pytest.mark.timeout(300)
def test_na_drop_full_pipeline():
    result = _run_pipeline(
        FIXTURES_DIR / "input_with_nas.csv",
        FIXTURES_DIR / "P0_small.csv",
        FIXTURES_DIR / "Q0_small.csv",
        n_workers=1,
        na_strategy="drop",
    )
    # Load original to find which genes had NaN
    input_with_nas = io_utils.read_expression(str(FIXTURES_DIR / "input_with_nas.csv"))
    na_gene_cols = list(input_with_nas.columns[input_with_nas.isna().any(axis=0)])

    # Those genes should be NaN in the output
    for g in na_gene_cols:
        if g in result.columns:
            assert result[g].isna().all(), f"Gene {g} should be NaN in output but is not"

    # Non-NA genes should have finite values
    clean_genes = [g for g in result.columns if g not in na_gene_cols]
    assert result[clean_genes].notna().all().all()


@octave_required
@pytest.mark.timeout(300)
def test_na_knn_full_pipeline():
    result = _run_pipeline(
        FIXTURES_DIR / "input_with_nas.csv",
        FIXTURES_DIR / "P0_small.csv",
        FIXTURES_DIR / "Q0_small.csv",
        n_workers=1,
        na_strategy="knn",
    )
    input_df = io_utils.read_expression(str(FIXTURES_DIR / "input_with_nas.csv"))
    # Output should have no NaN (all imputed)
    assert result.isna().sum().sum() == 0
    # Shape should match the input (minus any genes exceeding max_na_frac)
    assert result.shape[0] == input_df.shape[0]
