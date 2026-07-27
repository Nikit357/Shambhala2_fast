"""
Unit tests for the S0 harmonised output reuse path in run_shambhala_job.py.

Covers:
- _run_reuse_path: sample filtering, post-removal, r2 recording, upload, error handling
- Step 0 routing: Case A (all cached), Case B (reuse), Case B fallback to Case C
"""
import argparse
import pathlib
import sys
from unittest.mock import MagicMock, patch, call

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from harmonization_scripts.run_shambhala_job import _run_reuse_path
from harmonization_scripts.shambhala_bench_shared import (
    ALL_STRATEGIES,
    BATCH_COL,
    s3_key_exp,
    s3_key_prepared_ann,
)

# ── Shared helpers ─────────────────────────────────────────────────────────────

_N_SAMPLES = 20
_N_GENES = 10
_SAMPLE_IDS = [f"S{i}" for i in range(_N_SAMPLES)]
_GENE_IDS = [f"G{j}" for j in range(_N_GENES)]
_IMP = "strict"
_METHOD = "shambhala_P0std_Q0std"


def _make_exp(sample_ids: list[str], seed: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    return pd.DataFrame(
        rng.random((len(sample_ids), _N_GENES)),
        index=sample_ids,
        columns=_GENE_IDS,
    )


def _make_ann(sample_ids: list[str], n_batches: int = 4) -> pd.DataFrame:
    return pd.DataFrame(
        {
            BATCH_COL: [f"B{i % n_batches}" for i in range(len(sample_ids))],
            "Diagnosis_cell_type_unified": [f"D{i % 3}" for i in range(len(sample_ids))],
        },
        index=sample_ids,
    )


def _make_args(**overrides: object) -> argparse.Namespace:
    defaults: dict[str, object] = {"memory_limit_gb": 1.0, "skip_if_exists": True}
    defaults.update(overrides)
    return argparse.Namespace(**defaults)


def _build_output_specs(
    missing_strats: list[str],
    missing_post_rms: list[bool],
) -> tuple[list[tuple[str, bool, str]], list[tuple[str, bool, str]]]:
    all_specs: list[tuple[str, bool, str]] = [
        (strat, post_rm, s3_key_exp(strat, _IMP, _METHOD, post_rm))
        for strat in ALL_STRATEGIES
        for post_rm in [False, True]
    ]
    missing: list[tuple[str, bool, str]] = [
        (strat, post_rm, key)
        for strat, post_rm, key in all_specs
        if strat in missing_strats and post_rm in missing_post_rms
    ]
    return all_specs, missing


# ── _run_reuse_path unit tests ─────────────────────────────────────────────────

class TestRunReusePath:

    def test_derives_correct_subset_of_samples(self) -> None:
        exp_harmonized = _make_exp(_SAMPLE_IDS)
        ann_a = _make_ann(_SAMPLE_IDS[:10])
        all_specs, missing = _build_output_specs(["A_confirmed_bad"], [False])

        s3 = MagicMock()
        uploaded: dict[str, pd.DataFrame] = {}

        with (
            patch("harmonization_scripts.run_shambhala_job.download_exp_from_s3", return_value=exp_harmonized),
            patch("harmonization_scripts.run_shambhala_job.download_ann_from_s3", return_value=ann_a),
            patch("harmonization_scripts.run_shambhala_job.upload_exp_to_s3", side_effect=lambda df, s3c, k: uploaded.update({k: df.copy()})),
            patch("harmonization_scripts.run_shambhala_job.r2_batch", return_value=0.1),
            patch("harmonization_scripts.run_shambhala_job.free_gb", return_value=100.0),
        ):
            _run_reuse_path(
                s3=s3, imp=_IMP, method=_METHOD,
                missing_specs=missing, all_output_specs=all_specs,
                args=_make_args(), label="test",
            )

        assert len(uploaded) == 1
        key = missing[0][2]
        pd.testing.assert_frame_equal(uploaded[key], exp_harmonized.loc[_SAMPLE_IDS[:10]])

    def test_post_rm_removes_outlier_batch(self) -> None:
        exp_harmonized = _make_exp(_SAMPLE_IDS)
        ann_all = _make_ann(_SAMPLE_IDS)
        all_specs, missing = _build_output_specs(["A_confirmed_bad"], [True])

        s3 = MagicMock()
        uploaded: dict[str, pd.DataFrame] = {}
        outlier_batch = "B0"
        outlier_samples = [
            sid for sid, b in zip(_SAMPLE_IDS, ann_all[BATCH_COL]) if b == outlier_batch
        ]

        with (
            patch("harmonization_scripts.run_shambhala_job.download_exp_from_s3", return_value=exp_harmonized),
            patch("harmonization_scripts.run_shambhala_job.download_ann_from_s3", return_value=ann_all),
            patch("harmonization_scripts.run_shambhala_job.upload_exp_to_s3", side_effect=lambda df, s3c, k: uploaded.update({k: df.copy()})),
            patch("harmonization_scripts.run_shambhala_job.r2_batch", return_value=0.05),
            patch("harmonization_scripts.run_shambhala_job.free_gb", return_value=100.0),
            patch("harmonization_scripts.run_shambhala_job.identify_outlier_batches", return_value=[outlier_batch]),
        ):
            _run_reuse_path(
                s3=s3, imp=_IMP, method=_METHOD,
                missing_specs=missing, all_output_specs=all_specs,
                args=_make_args(), label="test",
            )

        key = missing[0][2]
        assert key in uploaded
        for sid in outlier_samples:
            assert sid not in uploaded[key].index

    def test_no_post_rm_when_flag_false(self) -> None:
        exp_harmonized = _make_exp(_SAMPLE_IDS)
        ann_all = _make_ann(_SAMPLE_IDS)
        all_specs, missing = _build_output_specs(["A_confirmed_bad"], [False])

        s3 = MagicMock()

        with (
            patch("harmonization_scripts.run_shambhala_job.download_exp_from_s3", return_value=exp_harmonized),
            patch("harmonization_scripts.run_shambhala_job.download_ann_from_s3", return_value=ann_all),
            patch("harmonization_scripts.run_shambhala_job.upload_exp_to_s3"),
            patch("harmonization_scripts.run_shambhala_job.r2_batch", return_value=0.1),
            patch("harmonization_scripts.run_shambhala_job.free_gb", return_value=100.0),
            patch("harmonization_scripts.run_shambhala_job.identify_outlier_batches") as mock_outlier,
        ):
            _run_reuse_path(
                s3=s3, imp=_IMP, method=_METHOD,
                missing_specs=missing, all_output_specs=all_specs,
                args=_make_args(), label="test",
            )

        mock_outlier.assert_not_called()

    def test_returns_all_24_rows(self) -> None:
        exp_harmonized = _make_exp(_SAMPLE_IDS)
        ann_a = _make_ann(_SAMPLE_IDS)
        all_specs, missing = _build_output_specs(["A_confirmed_bad"], [False])

        s3 = MagicMock()

        with (
            patch("harmonization_scripts.run_shambhala_job.download_exp_from_s3", return_value=exp_harmonized),
            patch("harmonization_scripts.run_shambhala_job.download_ann_from_s3", return_value=ann_a),
            patch("harmonization_scripts.run_shambhala_job.upload_exp_to_s3"),
            patch("harmonization_scripts.run_shambhala_job.r2_batch", return_value=0.1),
            patch("harmonization_scripts.run_shambhala_job.free_gb", return_value=100.0),
        ):
            rows = _run_reuse_path(
                s3=s3, imp=_IMP, method=_METHOD,
                missing_specs=missing, all_output_specs=all_specs,
                args=_make_args(), label="test",
            )

        assert len(rows) == 24

    def test_missing_rows_have_status_ok(self) -> None:
        exp_harmonized = _make_exp(_SAMPLE_IDS)
        ann_a = _make_ann(_SAMPLE_IDS)
        all_specs, missing = _build_output_specs(["A_confirmed_bad"], [False])

        s3 = MagicMock()

        with (
            patch("harmonization_scripts.run_shambhala_job.download_exp_from_s3", return_value=exp_harmonized),
            patch("harmonization_scripts.run_shambhala_job.download_ann_from_s3", return_value=ann_a),
            patch("harmonization_scripts.run_shambhala_job.upload_exp_to_s3"),
            patch("harmonization_scripts.run_shambhala_job.r2_batch", return_value=0.1),
            patch("harmonization_scripts.run_shambhala_job.free_gb", return_value=100.0),
        ):
            rows = _run_reuse_path(
                s3=s3, imp=_IMP, method=_METHOD,
                missing_specs=missing, all_output_specs=all_specs,
                args=_make_args(), label="test",
            )

        ok_rows = [r for r in rows if r["status"] == "ok"]
        cached_rows = [r for r in rows if r["status"] == "cached"]
        assert len(ok_rows) == 1
        assert ok_rows[0]["strat"] == "A_confirmed_bad"
        assert ok_rows[0]["post_rm"] is False
        assert len(cached_rows) == 23

    def test_cached_rows_have_nan_metrics(self) -> None:
        exp_harmonized = _make_exp(_SAMPLE_IDS)
        ann_a = _make_ann(_SAMPLE_IDS)
        all_specs, missing = _build_output_specs(["A_confirmed_bad"], [False])

        s3 = MagicMock()

        with (
            patch("harmonization_scripts.run_shambhala_job.download_exp_from_s3", return_value=exp_harmonized),
            patch("harmonization_scripts.run_shambhala_job.download_ann_from_s3", return_value=ann_a),
            patch("harmonization_scripts.run_shambhala_job.upload_exp_to_s3"),
            patch("harmonization_scripts.run_shambhala_job.r2_batch", return_value=0.1),
            patch("harmonization_scripts.run_shambhala_job.free_gb", return_value=100.0),
        ):
            rows = _run_reuse_path(
                s3=s3, imp=_IMP, method=_METHOD,
                missing_specs=missing, all_output_specs=all_specs,
                args=_make_args(), label="test",
            )

        for row in rows:
            if row["status"] == "cached":
                assert np.isnan(float(str(row["r2_batch"])))
                assert np.isnan(float(str(row["n_samples"])))

    def test_upload_failure_recorded_as_upload_failed(self) -> None:
        exp_harmonized = _make_exp(_SAMPLE_IDS)
        ann_a = _make_ann(_SAMPLE_IDS)
        all_specs, missing = _build_output_specs(["A_confirmed_bad"], [False])

        s3 = MagicMock()

        with (
            patch("harmonization_scripts.run_shambhala_job.download_exp_from_s3", return_value=exp_harmonized),
            patch("harmonization_scripts.run_shambhala_job.download_ann_from_s3", return_value=ann_a),
            patch("harmonization_scripts.run_shambhala_job.upload_exp_to_s3", side_effect=RuntimeError("S3 error")),
            patch("harmonization_scripts.run_shambhala_job.r2_batch", return_value=0.1),
            patch("harmonization_scripts.run_shambhala_job.free_gb", return_value=100.0),
        ):
            rows = _run_reuse_path(
                s3=s3, imp=_IMP, method=_METHOD,
                missing_specs=missing, all_output_specs=all_specs,
                args=_make_args(), label="test",
            )

        failed = [r for r in rows if r["status"] == "upload_failed"]
        assert len(failed) == 1
        assert failed[0]["strat"] == "A_confirmed_bad"

    def test_r2_failure_recorded_as_nan(self) -> None:
        exp_harmonized = _make_exp(_SAMPLE_IDS)
        ann_a = _make_ann(_SAMPLE_IDS)
        all_specs, missing = _build_output_specs(["A_confirmed_bad"], [False])

        s3 = MagicMock()

        with (
            patch("harmonization_scripts.run_shambhala_job.download_exp_from_s3", return_value=exp_harmonized),
            patch("harmonization_scripts.run_shambhala_job.download_ann_from_s3", return_value=ann_a),
            patch("harmonization_scripts.run_shambhala_job.upload_exp_to_s3"),
            patch("harmonization_scripts.run_shambhala_job.r2_batch", side_effect=RuntimeError("PCA error")),
            patch("harmonization_scripts.run_shambhala_job.free_gb", return_value=100.0),
        ):
            rows = _run_reuse_path(
                s3=s3, imp=_IMP, method=_METHOD,
                missing_specs=missing, all_output_specs=all_specs,
                args=_make_args(), label="test",
            )

        ok_rows = [r for r in rows if r["status"] == "ok"]
        assert len(ok_rows) == 1
        assert np.isnan(float(str(ok_rows[0]["r2_batch"])))
        assert np.isnan(float(str(ok_rows[0]["r2_diag"])))

    def test_raises_on_s0_download_failure(self) -> None:
        all_specs, missing = _build_output_specs(["A_confirmed_bad"], [False])
        s3 = MagicMock()

        with (
            patch("harmonization_scripts.run_shambhala_job.download_exp_from_s3", side_effect=RuntimeError("S3 unavailable")),
            patch("harmonization_scripts.run_shambhala_job.free_gb", return_value=100.0),
        ):
            with pytest.raises(RuntimeError, match="S3 unavailable"):
                _run_reuse_path(
                    s3=s3, imp=_IMP, method=_METHOD,
                    missing_specs=missing, all_output_specs=all_specs,
                    args=_make_args(), label="test",
                )

    def test_raises_on_insufficient_ram(self) -> None:
        all_specs, missing = _build_output_specs(["A_confirmed_bad"], [False])
        s3 = MagicMock()

        with patch("harmonization_scripts.run_shambhala_job.free_gb", return_value=0.5):
            with pytest.raises(RuntimeError, match="Insufficient RAM"):
                _run_reuse_path(
                    s3=s3, imp=_IMP, method=_METHOD,
                    missing_specs=missing, all_output_specs=all_specs,
                    args=_make_args(memory_limit_gb=1.0), label="test",
                )

    def test_downloads_annotations_only_for_missing_strategies(self) -> None:
        exp_harmonized = _make_exp(_SAMPLE_IDS)
        ann_a = _make_ann(_SAMPLE_IDS)
        all_specs, missing = _build_output_specs(["A_confirmed_bad"], [False])

        s3 = MagicMock()
        ann_download_keys: list[str] = []

        def record_ann_download(s3c: object, key: str) -> pd.DataFrame:
            ann_download_keys.append(key)
            return ann_a

        with (
            patch("harmonization_scripts.run_shambhala_job.download_exp_from_s3", return_value=exp_harmonized),
            patch("harmonization_scripts.run_shambhala_job.download_ann_from_s3", side_effect=record_ann_download),
            patch("harmonization_scripts.run_shambhala_job.upload_exp_to_s3"),
            patch("harmonization_scripts.run_shambhala_job.r2_batch", return_value=0.1),
            patch("harmonization_scripts.run_shambhala_job.free_gb", return_value=100.0),
        ):
            _run_reuse_path(
                s3=s3, imp=_IMP, method=_METHOD,
                missing_specs=missing, all_output_specs=all_specs,
                args=_make_args(), label="test",
            )

        expected_ann_key = s3_key_prepared_ann("A_confirmed_bad", _IMP)
        assert ann_download_keys == [expected_ann_key]


# ── Step 0 routing tests ───────────────────────────────────────────────────────

def _routing_argv(skip_if_exists: bool = True) -> list[str]:
    base = ["run_shambhala_job.py", "--imp", _IMP, "--method", _METHOD]
    if skip_if_exists:
        base.append("--skip-if-exists")
    return base


class TestStep0Routing:

    def test_case_a_exits_without_reuse_or_harmonization(self) -> None:
        all_keys = {
            s3_key_exp(strat, _IMP, _METHOD, post_rm)
            for strat in ALL_STRATEGIES
            for post_rm in [False, True]
        }

        with (
            patch("sys.argv", _routing_argv()),
            patch("boto3.client", return_value=MagicMock()),
            patch("harmonization_scripts.run_shambhala_job.s3_exists", side_effect=lambda s3c, k: k in all_keys),
            patch("harmonization_scripts.run_shambhala_job._write_sidecar"),
            patch("harmonization_scripts.run_shambhala_job._run_reuse_path") as mock_reuse,
        ):
            from harmonization_scripts import run_shambhala_job
            with pytest.raises(SystemExit) as exc:
                run_shambhala_job.main()

        assert exc.value.code == 0
        mock_reuse.assert_not_called()

    def test_case_b_calls_reuse_and_exits_zero(self) -> None:
        s0_key = s3_key_exp("S0_no_removal", _IMP, _METHOD, False)
        existing = {s0_key}

        fake_rows: list[dict[str, object]] = [{"strat": "S0_no_removal", "status": "ok"}]

        with (
            patch("sys.argv", _routing_argv()),
            patch("boto3.client", return_value=MagicMock()),
            patch("harmonization_scripts.run_shambhala_job.s3_exists", side_effect=lambda s3c, k: k in existing),
            patch("harmonization_scripts.run_shambhala_job._run_reuse_path", return_value=fake_rows) as mock_reuse,
            patch("harmonization_scripts.run_shambhala_job._write_sidecar"),
        ):
            from harmonization_scripts import run_shambhala_job
            with pytest.raises(SystemExit) as exc:
                run_shambhala_job.main()

        assert exc.value.code == 0
        mock_reuse.assert_called_once()

    def test_case_b_falls_through_to_case_c_on_reuse_failure(self) -> None:
        s0_key = s3_key_exp("S0_no_removal", _IMP, _METHOD, False)
        existing = {s0_key}

        with (
            patch("sys.argv", _routing_argv()),
            patch("boto3.client", return_value=MagicMock()),
            patch("harmonization_scripts.run_shambhala_job.s3_exists", side_effect=lambda s3c, k: k in existing),
            patch("harmonization_scripts.run_shambhala_job._run_reuse_path", side_effect=RuntimeError("download failed")),
            patch("harmonization_scripts.run_shambhala_job.free_gb", return_value=100.0),
            patch("harmonization_scripts.run_shambhala_job.s3_key_prepared_exp", return_value="fake/key"),
            patch("harmonization_scripts.run_shambhala_job.s3_key_prepared_ann", return_value="fake/ann"),
        ):
            from harmonization_scripts import run_shambhala_job
            with pytest.raises(SystemExit) as exc:
                run_shambhala_job.main()

        assert exc.value.code == 3

    def test_no_skip_if_exists_bypasses_step0(self) -> None:
        with (
            patch("sys.argv", _routing_argv(skip_if_exists=False)),
            patch("boto3.client", return_value=MagicMock()),
            patch("harmonization_scripts.run_shambhala_job.s3_exists", return_value=False),
            patch("harmonization_scripts.run_shambhala_job._run_reuse_path") as mock_reuse,
            patch("harmonization_scripts.run_shambhala_job.free_gb", return_value=100.0),
            patch("harmonization_scripts.run_shambhala_job.s3_key_prepared_exp", return_value="fake/key"),
            patch("harmonization_scripts.run_shambhala_job.s3_key_prepared_ann", return_value="fake/ann"),
        ):
            from harmonization_scripts import run_shambhala_job
            with pytest.raises(SystemExit) as exc:
                run_shambhala_job.main()

        mock_reuse.assert_not_called()
        assert exc.value.code == 3
