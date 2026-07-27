"""
Shared fixtures and helpers for speed-up integration tests.

All tests require Octave. They use the same 10-sample fixture as test_integration.py.
"""

import contextlib
import multiprocessing
import pathlib
import shutil
import sys
import time

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent.parent))

from shambhala import io_utils, na_handling, q_rescale, parallel as parallel_module


def assert_mean_rel_diff(
    metrics: dict[str, float],
    warn_threshold: float = 0.01,
    fail_threshold: float = 0.05,
    label: str = "",
) -> None:
    """
    Print all four metrics and assert two-tier mean_rel_diff thresholds.

    < warn_threshold  → STATUS: PASSED (production-safe)
    warn–fail range   → STATUS: WARNING (exploratory only, test still passes)
    > fail_threshold  → AssertionError (unacceptable distortion)
    """
    mrd = metrics["mean_rel_diff"]
    header = f"\n{label} metrics:" if label else "\nMetrics:"
    print(header)
    print(f"  max_abs_diff  = {metrics['max_abs_diff']:.4f}")
    print(f"  mean_abs_diff = {metrics['mean_abs_diff']:.4f}")
    print(f"  mean_rel_diff = {mrd:.4%}")
    print(f"  %within_1pct  = {metrics['pct_within_1pct']:.1f}%")

    if mrd < warn_threshold:
        print(f"  STATUS: PASSED — production-safe (< {warn_threshold:.0%})")
    elif mrd < fail_threshold:
        print(
            f"  STATUS: WARNING — exploratory only "
            f"({warn_threshold:.0%}–{fail_threshold:.0%} range). "
            f"Do not use for benchmarking."
        )

    assert mrd < fail_threshold, (
        f"mean_rel_diff={mrd:.4%} exceeds fail threshold {fail_threshold:.0%}"
    )

FIXTURES_DIR = pathlib.Path(__file__).parent.parent / "fixtures"
OCTAVE_DIR = pathlib.Path(__file__).resolve().parent.parent.parent / "octave"
OCTAVE_BIN = (
    "/opt/conda/bin/octave"
    if pathlib.Path("/opt/conda/bin/octave").exists()
    else (shutil.which("octave") or "octave")
)

_OCTAVE_AVAILABLE = (
    pathlib.Path("/opt/conda/bin/octave").exists() or shutil.which("octave") is not None
)
octave_required = pytest.mark.skipif(
    not _OCTAVE_AVAILABLE,
    reason="Octave not installed (expected at /opt/conda/bin/octave or on PATH)",
)


@pytest.fixture(autouse=True, scope="function")
def _kill_orphaned_workers():
    before = {p.pid for p in multiprocessing.active_children()}
    yield
    stragglers = [p for p in multiprocessing.active_children() if p.pid not in before]
    for p in stragglers:
        p.terminate()
    deadline = time.monotonic() + 3.0
    while stragglers and time.monotonic() < deadline:
        time.sleep(0.1)
        stragglers = [p for p in stragglers if p.is_alive()]
    for p in stragglers:
        p.kill()


def load_fixtures():
    """Return (input_df, p_df, q_df, expected_df) from the 10-sample fixture set."""
    input_df = io_utils.read_expression(str(FIXTURES_DIR / "input_10samples.csv"))
    p_df = io_utils.read_expression(str(FIXTURES_DIR / "P0_small.csv")).astype(float)
    q_df = io_utils.read_expression(str(FIXTURES_DIR / "Q0_small.csv")).astype(float)
    expected_df = io_utils.read_expression(str(FIXTURES_DIR / "expected_output.csv"))
    return input_df, p_df, q_df, expected_df


def run_pipeline(
    input_df: pd.DataFrame,
    p_df: pd.DataFrame,
    q_df: pd.DataFrame,
    n_workers: int = 1,
    na_strategy: str = "drop",
    random_seed: int = 42,
    precompute_qn_reference: bool = False,
    synthetic_cublock_p: bool = False,
    max_p_samples: int | None = None,
    precompute_cublock_clusters: bool = False,
    python_cublock: bool = False,
) -> pd.DataFrame:
    """
    Replicate run_shambhala.py main() logic with optional speed-up flags.

    Flags not yet implemented raise NotImplementedError; add implementation
    in each Phase as those approaches land.
    """
    pass  # python_cublock passed through to harmonize_parallel below

    input_T = input_df.T
    p_T = p_df.T
    q_T = q_df.T

    common_genes = input_T.index.intersection(p_T.index).intersection(q_T.index)
    input_T = input_T.loc[common_genes]
    p_T = p_T.loc[common_genes]
    q_T = q_T.loc[common_genes]

    # Approach C: centroid-based P subsampling
    if max_p_samples is not None and p_T.shape[1] > max_p_samples:
        p_df_aligned = p_T.T
        p_mean = p_df_aligned.mean(axis=0)
        distances = ((p_df_aligned - p_mean) ** 2).sum(axis=1)
        p_T = p_T[distances.nsmallest(max_p_samples).index]

    clean_df, na_mask = na_handling.apply_na_strategy(input_T.T, strategy=na_strategy)
    clean_T = clean_df.T
    p_T = p_T.loc[clean_T.index]
    q_df_filtered = q_T.loc[clean_T.index].T
    rm, rs = q_rescale.compute_q_statistics(q_df_filtered)

    p_for_octave = p_T

    # Approaches A + D: precompute QN reference from P
    if precompute_qn_reference or synthetic_cublock_p:
        p_sorted = np.sort(p_T.values, axis=0)
        p_qn_reference = p_sorted.mean(axis=1)

    if precompute_qn_reference:
        import qnorm as qnorm_lib  # noqa: PLC0415

        clean_T = qnorm_lib.quantile_normalize(clean_T, axis=1, target=p_qn_reference)

    if synthetic_cublock_p:
        p_for_octave = pd.DataFrame(
            p_qn_reference[:, np.newaxis],
            index=p_T.index,
            columns=["P_qn_reference"],
        )

    from shambhala.octave_bridge import precompute_clusters_from_p  # noqa: PLC0415

    fixed_clusters = None
    if precompute_cublock_clusters:
        fixed_clusters = precompute_clusters_from_p(
            p_T=p_for_octave,
            k=5,
            random_seed=random_seed,
        )

    harmonized_T = parallel_module.harmonize_parallel(
        input_df=clean_T,
        p_df=p_for_octave,
        rm=rm,
        rs=rs,
        k=5,
        n_workers=n_workers,
        octave_scripts_dir=str(OCTAVE_DIR),
        octave_bin=OCTAVE_BIN,
        timeout_s=300,
        random_seed=random_seed,
        fixed_clusters=fixed_clusters,
        python_cublock=python_cublock,
        disable_progress=True,
        skip_qn=precompute_qn_reference,
    )
    return na_handling.restore_na_genes(harmonized_T.T, na_mask)


@contextlib.contextmanager
def measure_wall_time():
    """Context manager that yields a list; after exit, list[0] is elapsed seconds."""
    result = []
    t0 = time.perf_counter()
    try:
        yield result
    finally:
        result.append(time.perf_counter() - t0)


def compare_to_expected(
    result: pd.DataFrame,
    expected: pd.DataFrame,
) -> dict[str, float]:
    """
    Compute faithfulness metrics between result and expected output.

    Returns dict with max_abs_diff, mean_abs_diff, mean_rel_diff, pct_within_1pct.
    """
    common_genes = result.columns.intersection(expected.columns)
    common_samples = result.index.intersection(expected.index)
    r = result.loc[common_samples, common_genes].values.astype(float)
    e = expected.loc[common_samples, common_genes].values.astype(float)

    diff = np.abs(r - e)
    denom = np.abs(e)
    rel_diff = np.where(denom > 0, diff / denom, 0.0)

    return {
        "max_abs_diff": float(diff.max()),
        "mean_abs_diff": float(diff.mean()),
        "mean_rel_diff": float(rel_diff.mean()),
        "pct_within_1pct": float((diff <= 0.01 * np.maximum(np.abs(e), 1e-9)).mean() * 100),
    }


def print_expression_comparison(
    result: pd.DataFrame,
    reference: pd.DataFrame,
    label: str = "",
    n_genes: int = 8,
    n_samples: int = 3,
) -> None:
    common_genes = result.columns.intersection(reference.columns)
    common_samples = result.index.intersection(reference.index)
    res = result.loc[common_samples, common_genes].astype(float)
    ref = reference.loc[common_samples, common_genes].astype(float)

    mean_abs_dev = (res - ref).abs().mean(axis=0)
    top_genes = mean_abs_dev.nlargest(n_genes).index.tolist()

    samples = sorted(common_samples.tolist())[:n_samples]

    header = f"--- Expression comparison: {label} (first {n_samples} samples, top {n_genes} changed genes) ---"
    print(f"\n{header}")

    col_w = 24
    gene_w = 16
    sample_header = "".join(f"{s:<{col_w}}" for s in samples)
    print(f"{'Gene':<{gene_w}}{sample_header}")

    subheader = "".join(f"{'ref':>10}  {'result':>10}  " for _ in samples)
    print(f"{'':<{gene_w}}{subheader}")

    for gene in top_genes:
        row = f"{gene:<{gene_w}}"
        for s in samples:
            rv = ref.at[s, gene]
            av = res.at[s, gene]
            row += f"{rv:>10.2f}  {av:>10.2f}  "
        print(row)
    print()
