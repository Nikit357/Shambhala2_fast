"""
Derive I_rare_batches_removed outputs from existing S0 harmonized datasets on S3.

Shambhala normalizes each sample independently. The strategy I outputs can be
derived by filtering rows from the existing S0 harmonized matrix without re-running
Shambhala.

Reads from S3:
    FL_batch_correction/prepared/S0_no_removal__strict__ann.tsv.gz
    FL_batch_correction/exp/S0_no_removal__{imp}__{method}__post0.tsv.gz  (per job)

Writes to S3:
    FL_batch_correction/prepared/I_rare_batches_removed__{imp}__ann.tsv.gz  (3 imps)
    FL_batch_correction/exp/I_rare_batches_removed__{imp}__{method}__post{0|1}.tsv.gz  (108)
    FL_batch_correction/shambhala_metrics/{imp}__{method}.json  (updated, 54)
    FL_batch_correction/shambhala_metrics.csv  (re-consolidated)

Exit codes
----------
0  Completed (some jobs may be skipped if S0 post0 is absent).
1  S0 annotation not found — cannot compute strategy I sample IDs.

Usage
-----
python harmonization_scripts/derive_rare_batches_outputs.py \\
    --skip-if-exists \\
    --log-level INFO
"""
from __future__ import annotations

import argparse
import gc
import io
import json
import logging
import sys
import time
import traceback
from pathlib import Path

import boto3
import pandas as pd

_THIS_DIR = Path(__file__).resolve().parent
_SHAMBHALA_ROOT = _THIS_DIR.parent
if str(_SHAMBHALA_ROOT) not in sys.path:
    sys.path.insert(0, str(_SHAMBHALA_ROOT))

from harmonization_scripts.shambhala_bench_shared import (
    ALL_IMPUTATION,
    ALL_METHODS,
    BATCH_COL,
    BIO_COL,
    S3_BUCKET,
    download_ann_from_s3,
    download_exp_from_s3,
    identify_outlier_batches,
    r2_batch,
    s3_exists,
    s3_key_exp,
    s3_key_prepared_ann,
    s3_key_sidecar,
    s3_key_shambhala_metrics_csv,
    upload_exp_to_s3,
)

STRAT_I = "I_rare_batches_removed"
MIN_BATCH_SIZE = 50


def _ts() -> str:
    return time.strftime("%H:%M:%S")


def _log(msg: str) -> None:
    print(f"[{_ts()}] {msg}", flush=True)


def _compute_strat_i_ids(ann_s0: pd.DataFrame) -> pd.Index:
    counts = ann_s0[BATCH_COL].value_counts()
    large = counts[counts >= MIN_BATCH_SIZE].index
    return ann_s0[ann_s0[BATCH_COL].isin(large)].index


def _read_sidecar(s3_client: object, imp: str, method: str) -> list[dict]:
    key = s3_key_sidecar(imp, method)
    try:
        body = s3_client.get_object(Bucket=S3_BUCKET, Key=key)["Body"].read()  # type: ignore[attr-defined]
        return json.loads(body)
    except Exception:
        return []


def _write_sidecar(s3_client: object, imp: str, method: str, rows: list[dict]) -> None:
    key = s3_key_sidecar(imp, method)
    buf = io.BytesIO(json.dumps(rows, indent=2).encode())
    try:
        s3_client.upload_fileobj(buf, S3_BUCKET, key)  # type: ignore[attr-defined]
    except Exception:
        _log(f"WARNING: could not upload sidecar {key}")


def _reconsolidate_metrics(s3_client: object) -> None:
    rows: list[dict] = []
    for imp in ALL_IMPUTATION:
        for method in ALL_METHODS:
            key = s3_key_sidecar(imp, method)
            try:
                body = s3_client.get_object(Bucket=S3_BUCKET, Key=key)["Body"].read()  # type: ignore[attr-defined]
                rows.extend(json.loads(body))
            except Exception:
                pass
    if not rows:
        _log("No sidecar rows found — shambhala_metrics.csv not updated.")
        return
    df = pd.DataFrame(rows).drop_duplicates(
        subset=["strat", "imp", "method", "post_rm"], keep="last"
    )
    key = s3_key_shambhala_metrics_csv()
    s3_client.upload_fileobj(io.BytesIO(df.to_csv(index=False).encode()), S3_BUCKET, key)  # type: ignore[attr-defined]
    _log(f"Consolidated metrics uploaded: s3://{S3_BUCKET}/{key} ({len(df)} rows)")


def main() -> None:
    logging.basicConfig(level=logging.INFO, stream=sys.stdout)
    parser = argparse.ArgumentParser(
        description=f"Derive {STRAT_I} S3 outputs from existing S0 harmonized datasets."
    )
    parser.add_argument("--skip-if-exists", action="store_true", dest="skip_if_exists")
    parser.add_argument("--imps", default=None, help="Comma-separated imputation subset (default: all 3).")
    parser.add_argument("--methods", default=None, help="Comma-separated method subset (default: all 18).")
    parser.add_argument(
        "--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING"], dest="log_level"
    )
    args = parser.parse_args()

    logging.getLogger().setLevel(getattr(logging, args.log_level))

    imps: list[str] = args.imps.split(",") if args.imps else ALL_IMPUTATION
    methods: list[str] = args.methods.split(",") if args.methods else ALL_METHODS

    s3 = boto3.client("s3")

    # ── Step 1: Compute strategy I sample IDs from S0 annotation ─────────────────
    ann_s0_key = s3_key_prepared_ann("S0_no_removal", "strict")
    if not s3_exists(s3, ann_s0_key):
        print(f"S0 annotation not found: {ann_s0_key}", file=sys.stderr)
        sys.exit(1)
    _log("Downloading S0 annotation ...")
    ann_s0 = download_ann_from_s3(s3, ann_s0_key)
    strat_i_ids = _compute_strat_i_ids(ann_s0)
    ann_i = ann_s0.loc[strat_i_ids]
    _log(
        f"Strategy I: {len(strat_i_ids)} / {len(ann_s0)} samples survive "
        f"MIN_BATCH_SIZE={MIN_BATCH_SIZE} filter."
    )
    del ann_s0
    gc.collect()

    # ── Step 2: Upload strategy I annotation for all imputation keys ──────────────
    for imp in ALL_IMPUTATION:
        key = s3_key_prepared_ann(STRAT_I, imp)
        if args.skip_if_exists and s3_exists(s3, key):
            _log(f"Annotation already on S3 — skipping: {STRAT_I}__{imp}")
            continue
        _log(f"Uploading annotation: {STRAT_I}__{imp} ...")
        upload_exp_to_s3(ann_i, s3, key)

    # ── Step 3: Derive per-(imp, method) outputs ──────────────────────────────────
    n_total = len(imps) * len(methods)
    n_done = n_derived = n_skipped = n_no_s0 = 0

    for imp in imps:
        for method in methods:
            n_done += 1
            tag = f"{imp}×{method}"

            if args.skip_if_exists:
                key0 = s3_key_exp(STRAT_I, imp, method, False)
                key1 = s3_key_exp(STRAT_I, imp, method, True)
                if s3_exists(s3, key0) and s3_exists(s3, key1):
                    _log(f"[{n_done}/{n_total}] {tag}: already derived — skipping.")
                    n_skipped += 1
                    continue

            s0_key = s3_key_exp("S0_no_removal", imp, method, False)
            if not s3_exists(s3, s0_key):
                _log(f"[{n_done}/{n_total}] {tag}: S0 post0 absent — skipping.")
                n_no_s0 += 1
                continue

            _log(f"[{n_done}/{n_total}] {tag}: downloading S0 post0 ...")
            try:
                exp_s0 = download_exp_from_s3(s3, s0_key)
            except Exception:
                _log(f"  Download failed:\n{traceback.format_exc()}")
                continue

            common_ids = exp_s0.index.intersection(strat_i_ids)
            exp_i = exp_s0.loc[common_ids]
            ann_i_job = ann_i.loc[ann_i.index.intersection(common_ids)]
            _log(f"  {len(common_ids)} samples after filter ({exp_i.shape[1]} genes)")
            del exp_s0
            gc.collect()

            job_rows: list[dict] = []
            for post_rm in [False, True]:
                if post_rm:
                    exp_out = exp_i.copy()
                    ann_out = ann_i_job.copy()
                    try:
                        outliers = identify_outlier_batches(exp_out, ann_out)
                        if outliers:
                            _log(f"  post1: removing {len(outliers)} outlier batch(es): {outliers}")
                            ann_out = ann_out[~ann_out[BATCH_COL].isin(outliers)]
                            exp_out = exp_out.loc[ann_out.index]
                        else:
                            _log("  post1: no outlier batches found.")
                    except Exception:
                        _log(f"  post1 outlier detection failed:\n{traceback.format_exc()}")
                else:
                    exp_out = exp_i
                    ann_out = ann_i_job

                try:
                    r2_b = r2_batch(exp_out, ann_out, batch_col=BATCH_COL)
                    r2_d = r2_batch(exp_out, ann_out, batch_col=BIO_COL)
                except Exception:
                    r2_b = r2_d = float("nan")

                out_key = s3_key_exp(STRAT_I, imp, method, post_rm)
                try:
                    upload_exp_to_s3(exp_out, s3, out_key)
                    status = "derived"
                    _log(
                        f"  post_rm={post_rm}: uploaded ({len(ann_out)} samples, "
                        f"r2_batch={r2_b:.3f})"
                    )
                except Exception:
                    status = "upload_failed"
                    _log(f"  post_rm={post_rm}: upload FAILED:\n{traceback.format_exc()}")

                job_rows.append({
                    "strat": STRAT_I, "imp": imp, "method": method, "post_rm": post_rm,
                    "r2_batch": r2_b, "r2_diag": r2_d,
                    "n_samples": len(ann_out), "n_genes": exp_out.shape[1],
                    "status": status,
                })

                if post_rm:
                    del exp_out, ann_out

            del exp_i, ann_i_job
            gc.collect()
            n_derived += 1

            existing = _read_sidecar(s3, imp, method)
            existing = [r for r in existing if r.get("strat") != STRAT_I]
            _write_sidecar(s3, imp, method, existing + job_rows)

    _log(f"Derived: {n_derived}, skipped (already on S3): {n_skipped}, no S0 post0: {n_no_s0}.")

    if n_derived > 0 or n_skipped > 0:
        _log("Re-consolidating shambhala_metrics.csv ...")
        _reconsolidate_metrics(s3)


if __name__ == "__main__":
    main()
