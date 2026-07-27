"""
Unit tests for shambhala/octave_bridge.py.
"""

import io
import pathlib
import shutil
import sys
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import queue

from shambhala.octave_bridge import (
    OctaveExecutionError,
    OctaveOutputError,
    _format_pool_as_tsv,
    _is_progress_line,
    _parse_octave_stdout,
    _parse_progress_line,
    run_octave_normalize,
    run_octave_normalize_streaming,
)
from shambhala.progress_display import ProgressEvent

FIXTURES_DIR = pathlib.Path(__file__).parent / "fixtures"
OCTAVE_DIR = pathlib.Path(__file__).resolve().parent.parent / "octave"

# Detect octave binary — check /opt/conda/bin/octave first, then PATH
_OCTAVE_BIN = "/opt/conda/bin/octave" if pathlib.Path("/opt/conda/bin/octave").exists() else shutil.which("octave")
_OCTAVE_AVAILABLE = _OCTAVE_BIN is not None

octave_required = pytest.mark.skipif(
    not _OCTAVE_AVAILABLE,
    reason="Octave not installed (expected at /opt/conda/bin/octave or on PATH)",
)


def _make_pool_df(n_genes: int = 20, n_samples: int = 3, n_p: int = 5) -> pd.DataFrame:
    rng = np.random.default_rng(7)
    data = rng.random((n_genes, n_samples + n_p)) * 500
    genes = [f"G{i}" for i in range(n_genes)]
    samples = [f"S{i}" for i in range(n_samples)] + [f"P{i}" for i in range(n_p)]
    return pd.DataFrame(data, index=genes, columns=samples)


@octave_required
def test_bridge_returns_correct_shape():
    pool_df = _make_pool_df(n_genes=30, n_samples=5, n_p=5)
    result = run_octave_normalize(
        pool_df=pool_df,
        nh=5,
        np_=5,
        k=5,
        octave_scripts_dir=str(OCTAVE_DIR),
        octave_bin=_OCTAVE_BIN,
        timeout_s=120,
        random_seed=42,
    )
    assert result.shape == (30, 5)


@octave_required
def test_bridge_symbol_index_correct():
    pool_df = _make_pool_df(n_genes=15, n_samples=3, n_p=5)
    result = run_octave_normalize(
        pool_df=pool_df,
        nh=3,
        np_=5,
        k=5,
        octave_scripts_dir=str(OCTAVE_DIR),
        octave_bin=_OCTAVE_BIN,
        timeout_s=120,
        random_seed=42,
    )
    assert list(result.index) == list(pool_df.index)
    assert list(result.columns) == ["S0", "S1", "S2"]


def test_bridge_raises_on_nonzero_exit():
    pool_df = _make_pool_df()
    with pytest.raises(OctaveExecutionError):
        run_octave_normalize(
            pool_df=pool_df,
            nh=3,
            np_=5,
            k=5,
            octave_scripts_dir=str(OCTAVE_DIR),
            octave_bin="/nonexistent/octave",
            timeout_s=10,
        )


def test_bridge_error_message_contains_stderr():
    pool_df = _make_pool_df()
    fake_result = MagicMock()
    fake_result.returncode = 1
    fake_result.stderr = "error: undefined symbol NH"
    fake_result.stdout = ""
    with patch("subprocess.run", return_value=fake_result):
        with pytest.raises(OctaveExecutionError) as exc_info:
            run_octave_normalize(
                pool_df=pool_df,
                nh=3,
                np_=5,
                k=5,
                octave_scripts_dir=str(OCTAVE_DIR),
            )
    assert "error: undefined symbol NH" in str(exc_info.value)


def test_bridge_raises_on_unparseable_stdout():
    pool_df = _make_pool_df(n_genes=5, n_samples=3, n_p=5)
    fake_result = MagicMock()
    fake_result.returncode = 0
    fake_result.stdout = "garbage output that cannot be parsed into a matrix"
    fake_result.stderr = ""
    with patch("subprocess.run", return_value=fake_result):
        with pytest.raises(OctaveOutputError):
            run_octave_normalize(
                pool_df=pool_df,
                nh=3,
                np_=5,
                k=5,
                octave_scripts_dir=str(OCTAVE_DIR),
            )


def test_tsv_formatting_round_trip():
    pool_df = _make_pool_df(n_genes=5, n_samples=2, n_p=2)
    tsv_str = _format_pool_as_tsv(pool_df)
    parsed = pd.read_csv(io.StringIO(tsv_str), sep="\t", index_col=0)
    parsed.index.name = None  # TSV has "SYMBOL" as index column name; pool_df has None
    pd.testing.assert_frame_equal(parsed, pool_df, check_dtype=False, atol=1e-10)


# ── _is_progress_line / _parse_progress_line ──────────────────────────────────


def test_is_progress_line_true():
    assert _is_progress_line("Harmonizing sample 3 out of 10") is True


def test_is_progress_line_true_extra_spaces():
    assert _is_progress_line("Harmonizing sample  3  out of  10") is True


def test_is_progress_line_false_data_line():
    assert _is_progress_line("TP53 1.234 5.678") is False


def test_is_progress_line_false_empty():
    assert _is_progress_line("") is False


def test_is_progress_line_false_partial():
    assert _is_progress_line("Harmonizing") is False


def test_parse_progress_line_valid():
    x, y = _parse_progress_line("Harmonizing sample 3 out of 10")
    assert x == 3
    assert y == 10


def test_parse_progress_line_extra_spaces():
    x, y = _parse_progress_line("Harmonizing sample  7  out of  20")
    assert x == 7
    assert y == 20


def test_parse_progress_line_invalid_raises():
    with pytest.raises(ValueError, match="Not a progress line"):
        _parse_progress_line("TP53 1.234 5.678")


# ── run_octave_normalize_streaming ────────────────────────────────────────────


def test_streaming_bad_binary_raises():
    pool_df = _make_pool_df()
    with pytest.raises(OctaveExecutionError, match="not found"):
        run_octave_normalize_streaming(
            pool_df=pool_df,
            nh=3,
            np_=5,
            k=5,
            octave_scripts_dir=str(OCTAVE_DIR),
            octave_bin="nonexistent_binary_xyz",
            timeout_s=10,
        )


@octave_required
def test_streaming_no_queue_matches_blocking():
    pool_df = _make_pool_df(n_genes=30, n_samples=5, n_p=5)
    result_blocking = run_octave_normalize(
        pool_df=pool_df, nh=5, np_=5, k=5,
        octave_scripts_dir=str(OCTAVE_DIR), octave_bin=_OCTAVE_BIN,
        timeout_s=120, random_seed=42,
    )
    result_streaming = run_octave_normalize_streaming(
        pool_df=pool_df, nh=5, np_=5, k=5,
        octave_scripts_dir=str(OCTAVE_DIR), octave_bin=_OCTAVE_BIN,
        timeout_s=120, random_seed=42,
        progress_queue=None,
    )
    np.testing.assert_allclose(result_blocking.values, result_streaming.values, atol=1e-6)


@octave_required
def test_streaming_sends_progress_events():
    pool_df = _make_pool_df(n_genes=20, n_samples=3, n_p=5)
    q: queue.Queue[ProgressEvent | None] = queue.Queue()
    run_octave_normalize_streaming(
        pool_df=pool_df, nh=3, np_=5, k=5,
        octave_scripts_dir=str(OCTAVE_DIR), octave_bin=_OCTAVE_BIN,
        timeout_s=120, random_seed=42,
        progress_queue=q, worker_idx=2,
    )
    events: list[ProgressEvent] = []
    while not q.empty():
        item = q.get_nowait()
        if item is not None:
            events.append(item)

    # Expect 3 non-final + 1 final = 4 events total
    non_final = [e for e in events if not e.is_final]
    final = [e for e in events if e.is_final]
    assert len(non_final) == 3
    assert len(final) == 1
    assert final[0].worker_idx == 2
    assert final[0].sample_done == 3
    assert final[0].sample_total == 3
    # sample_done values should be 1, 2, 3
    assert [e.sample_done for e in non_final] == [1, 2, 3]


def test_random_seed_injected_in_preamble():
    pool_df = _make_pool_df(n_genes=5, n_samples=2, n_p=2)
    fake_result = MagicMock()
    fake_result.returncode = 0
    # Return a valid stdout with correct shape
    genes = list(pool_df.index)
    fake_stdout_lines = [f"{g} 1.0 2.0" for g in genes]
    fake_result.stdout = "\n".join(fake_stdout_lines)
    fake_result.stderr = ""

    captured_cmd = []
    def fake_run(cmd, **kwargs):
        captured_cmd.extend(cmd)
        return fake_result

    with patch("subprocess.run", side_effect=fake_run):
        run_octave_normalize(
            pool_df=pool_df,
            nh=2,
            np_=2,
            k=5,
            octave_scripts_dir=str(OCTAVE_DIR),
            random_seed=42,
        )

    eval_str = " ".join(captured_cmd)
    assert "rand('state', 42)" in eval_str
    assert "randn('state', 42)" in eval_str
