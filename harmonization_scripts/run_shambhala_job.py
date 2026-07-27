"""
Shambhala cross-product worker: runs Shambhala on the full S0 dataset for one
(imp × method) pair, then writes all 28 (strategy × post_rm) outputs to S3.

Three execution paths (selected automatically when --skip-if-exists is set):
  Case A — all 28 outputs already on S3: early exit, no downloads.
  Case B — S0 harmonized (S0_no_removal__post0) present, some outputs missing:
           download S0 file, filter rows per strategy, upload missing outputs.
           No Octave invocation.
  Case C — S0 harmonized missing (or --skip-if-exists absent): full harmonization
           via Octave + download of S0 prepared data and calibration P/Q.

Reads from S3:
    FL_batch_correction/prepared/S0_no_removal__{imp}__exp.tsv.gz           (Case C)
    FL_batch_correction/prepared/S0_no_removal__{imp}__ann.tsv.gz           (Case C)
    FL_batch_correction/prepared/{strat}__{imp}__ann.tsv.gz   (×14 strategies)
    FL_batch_correction/calibration_datasets/{P_FILE}.csv                   (Case C)
    FL_batch_correction/calibration_datasets/{Q_FILE}.csv                   (Case C)
    FL_batch_correction/exp/S0_no_removal__{imp}__{method}__post0.tsv.gz    (Case B)

Writes to S3:
    FL_batch_correction/exp/{strat}__{imp}__{method}__post{0|1}.tsv.gz  (×28)
    FL_batch_correction/shambhala_metrics/{imp}__{method}.json

Exit codes
----------
0  Success (or all outputs already cached).
1  Shambhala harmonization failed.
2  Insufficient RAM (pre-flight check).
3  Prepared S0 dataset not found in S3.

Usage
-----
python run_shambhala_job.py \\
    --imp knn \\
    --method shambhala_P0std_Q0std \\
    [--skip-if-exists] \\
    [--out-json /tmp/shambhala_jobs/knn__shambhala_P0std_Q0std.json] \\
    [--memory-limit-gb 8.0] \\
    [--n-shambhala-workers 5] \\
    [--octave-bin octave] \\
    [--random-seed 42] \\
    [--shambhala-timeout-s 7200] \\
    [--na-strategy drop]
"""
from __future__ import annotations

import argparse
import gc
import json
import logging
import pathlib
import sys
import time
import traceback

import boto3
import numpy as np
import pandas as pd

# Allow running from the harmonization_scripts/ directory directly
_THIS_DIR = pathlib.Path(__file__).resolve().parent
_SHAMBHALA_ROOT = _THIS_DIR.parent
if str(_SHAMBHALA_ROOT) not in sys.path:
    sys.path.insert(0, str(_SHAMBHALA_ROOT))

from harmonization_scripts.shambhala_bench_shared import (
    ALL_STRATEGIES,
    BATCH_COL,
    BIO_COL,
    P_DATASETS,
    Q_DATASETS,
    SHAMBHALA_VARIANTS,
    download_ann_from_s3,
    download_calib_df,
    download_exp_from_s3,
    free_gb,
    identify_outlier_batches,
    r2_batch,
    s3_exists,
    s3_key_exp,
    s3_key_prepared_ann,
    s3_key_prepared_exp,
    s3_key_sidecar,
    upload_exp_to_s3,
    S3_BUCKET,
)
from shambhala import na_handling, q_rescale, parallel as parallel_module

_OCTAVE_DIR = _SHAMBHALA_ROOT / "octave"


def _ts() -> str:
    return time.strftime("%H:%M:%S")


def _log(label: str, msg: str) -> None:
    print(f"[{_ts()}][{label}] {msg}", flush=True)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Shambhala cross-product worker: one (imp × method) → 28 S3 outputs."
    )
    parser.add_argument("--imp", required=True, help="Imputation method (strict/knn/softimpute).")
    parser.add_argument("--method", required=True, help="Shambhala variant key (e.g. shambhala_P0std_Q0std).")
    parser.add_argument("--skip-if-exists", action="store_true", dest="skip_if_exists")
    parser.add_argument("--out-json", default=None, dest="out_json", metavar="PATH")
    parser.add_argument("--memory-limit-gb", type=float, default=8.0, dest="memory_limit_gb")
    parser.add_argument("--n-shambhala-workers", type=int, default=5, dest="n_shambhala_workers")
    parser.add_argument("--octave-bin", default="octave", dest="octave_bin")
    parser.add_argument("--random-seed", type=int, default=None, dest="random_seed")
    parser.add_argument("--shambhala-timeout-s", type=int, default=7200, dest="shambhala_timeout_s")
    parser.add_argument("--na-strategy", choices=["drop", "knn"], default="drop", dest="na_strategy")
    parser.add_argument(
        "--precompute-qn-reference", action="store_true", default=False,
        dest="precompute_qn_reference",
        help="Precompute QN reference from P once; apply Python-side QN per sample. ~10-20× speedup.",
    )
    parser.add_argument(
        "--synthetic-cublock-p", action="store_true", default=False,
        dest="synthetic_cublock_p",
        help="Pass 1-column synthetic P to Octave CuBlock instead of full P.",
    )
    parser.add_argument(
        "--max-p-samples", type=int, default=None, dest="max_p_samples", metavar="N",
        help="Subsample P to N centroid-closest samples (default: disabled).",
    )
    parser.add_argument(
        "--precompute-cublock-clusters", action="store_true", default=False,
        dest="precompute_cublock_clusters",
        help="Precompute k-means gene clusters from P once; reuse per sample (not yet implemented).",
    )
    parser.add_argument(
        "--python-cublock", action="store_true", default=False,
        dest="python_cublock",
        help="Use Python CuBlock instead of Octave (not yet implemented).",
    )
    return parser


def _make_failed_rows(imp: str, method: str, status: str = "failed") -> list[dict]:
    return [
        {
            "strat": strat,
            "imp": imp,
            "method": method,
            "post_rm": post_rm,
            "r2_batch": float("nan"),
            "r2_diag": float("nan"),
            "n_samples": float("nan"),
            "n_genes": float("nan"),
            "status": status,
        }
        for strat in ALL_STRATEGIES
        for post_rm in [False, True]
    ]


def _write_sidecar(
    rows: list[dict],
    s3_client: object,
    imp: str,
    method: str,
    out_json: str | None,
) -> None:
    import io as _io
    key = s3_key_sidecar(imp, method)
    body = json.dumps(rows, indent=2).encode()
    buf = _io.BytesIO(body)
    try:
        s3_client.upload_fileobj(buf, S3_BUCKET, key)  # type: ignore[attr-defined]
    except Exception:
        print(f"[{_ts()}] WARNING: could not upload sidecar to S3: {key}", flush=True)
    if out_json:
        pathlib.Path(out_json).parent.mkdir(parents=True, exist_ok=True)
        with open(out_json, "w") as fh:
            json.dump(rows, fh, indent=2)


def _run_reuse_path(
    s3: object,
    imp: str,
    method: str,
    missing_specs: list[tuple[str, bool, str]],
    all_output_specs: list[tuple[str, bool, str]],
    args: argparse.Namespace,
    label: str,
) -> list[dict[str, object]]:
    """
    Derive missing (strat × post_rm) outputs from the pre-existing S0 harmonized
    matrix on S3. Downloads the S0 output, filters rows per strategy, uploads
    missing files. No Octave invocation.

    Raises on S0 download failure so the caller can fall back to full harmonization.
    """
    ram = free_gb()
    if ram < args.memory_limit_gb:
        raise RuntimeError(
            f"Insufficient RAM for reuse path: {ram:.1f} GB free < "
            f"{args.memory_limit_gb:.1f} GB required."
        )
    _log(label, f"Reuse path RAM pre-flight OK: {ram:.1f} GB free.")

    t_total = time.time()

    s0_key = s3_key_exp("S0_no_removal", imp, method, False)
    _log(label, f"Downloading S0 harmonized matrix from {s0_key} ...")
    t0 = time.time()
    exp_harmonized = download_exp_from_s3(s3, s0_key)
    _log(
        label,
        f"S0 harmonized: {exp_harmonized.shape[0]} samples × "
        f"{exp_harmonized.shape[1]} genes ({time.time()-t0:.0f}s)",
    )

    strats_needed = {strat for strat, _, _ in missing_specs}
    _log(label, f"Downloading annotations for {len(strats_needed)} strategy/ies ...")
    t0 = time.time()
    ann_strat: dict[str, pd.DataFrame] = {}
    for strat in strats_needed:
        ann_strat[strat] = download_ann_from_s3(s3, s3_key_prepared_ann(strat, imp))
    _log(label, f"Annotations ready ({time.time()-t0:.0f}s).")

    missing_keys = {key for _, _, key in missing_specs}
    rows: list[dict[str, object]] = []
    n_done = 0

    for strat, post_rm, key in missing_specs:
        n_done += 1
        tag = f"{strat}×post_rm={post_rm}"

        strat_sample_ids = exp_harmonized.index.intersection(ann_strat[strat].index)
        exp_strat = exp_harmonized.loc[strat_sample_ids]
        ann_s = ann_strat[strat].loc[strat_sample_ids]
        _log(label, f"[{n_done}/{len(missing_specs)}] {tag}: {len(ann_s)} samples after filter")

        if post_rm:
            try:
                outliers = identify_outlier_batches(exp_strat, ann_s)
                if outliers:
                    _log(label, f"  Removing {len(outliers)} outlier batch(es): {outliers}")
                    ann_s = ann_s[~ann_s[BATCH_COL].isin(outliers)]
                    exp_strat = exp_strat.loc[ann_s.index]
                else:
                    _log(label, "  No outlier batches found.")
            except Exception:
                _log(label, f"  Post-removal failed (continuing):\n{traceback.format_exc()}")

        try:
            r2_b = r2_batch(exp_strat, ann_s, batch_col=BATCH_COL)
            r2_d = r2_batch(exp_strat, ann_s, batch_col=BIO_COL)
        except Exception:
            r2_b = r2_d = float("nan")

        _log(label, f"  r2_batch={r2_b:.3f} r2_diag={r2_d:.3f} — uploading {tag} ...")
        try:
            upload_exp_to_s3(exp_strat, s3, key)
            status = "ok"
            _log(label, "  Upload complete.")
        except Exception:
            status = "upload_failed"
            _log(label, f"  Upload FAILED (continuing):\n{traceback.format_exc()}")

        rows.append({
            "strat": strat,
            "imp": imp,
            "method": method,
            "post_rm": post_rm,
            "r2_batch": r2_b,
            "r2_diag": r2_d,
            "n_samples": len(ann_s),
            "n_genes": exp_strat.shape[1],
            "status": status,
        })
        del exp_strat, ann_s
        gc.collect()

    for strat, post_rm, key in all_output_specs:
        if key not in missing_keys:
            rows.append({
                "strat": strat,
                "imp": imp,
                "method": method,
                "post_rm": post_rm,
                "r2_batch": float("nan"),
                "r2_diag": float("nan"),
                "n_samples": float("nan"),
                "n_genes": float("nan"),
                "status": "cached",
            })

    elapsed = time.time() - t_total
    n_ok = sum(1 for r in rows if r["status"] == "ok")
    n_cached = sum(1 for r in rows if r["status"] == "cached")
    _log(
        label,
        f"Reuse path done: {n_ok} derived, {n_cached} already cached. "
        f"Total elapsed: {elapsed:.0f}s",
    )
    return rows


def main() -> None:
    logging.basicConfig(
        level=logging.INFO,
        stream=sys.stdout,
        format="%(asctime)s %(name)s %(levelname)s %(message)s",
    )
    args = _build_parser().parse_args()
    imp = args.imp
    method = args.method
    label = f"{imp}×{method}"

    if method not in SHAMBHALA_VARIANTS:
        print(
            f"[{_ts()}] Unknown method: {method!r}. Valid: {sorted(SHAMBHALA_VARIANTS)[:3]} ...",
            file=sys.stderr,
        )
        sys.exit(1)

    s3 = boto3.client("s3")

    # ── Step 0: Cache and reuse check ─────────────────────────────────────────
    if args.skip_if_exists:
        all_output_specs: list[tuple[str, bool, str]] = [
            (strat, post_rm, s3_key_exp(strat, imp, method, post_rm))
            for strat in ALL_STRATEGIES
            for post_rm in [False, True]
        ]
        missing_specs = [
            (strat, post_rm, key)
            for strat, post_rm, key in all_output_specs
            if not s3_exists(s3, key)
        ]

        if not missing_specs:
            # Case A: all 28 outputs already cached
            _log(label, "All 28 outputs already on S3 — cached.")
            rows = _make_failed_rows(imp, method, status="cached")
            _write_sidecar(rows, s3, imp, method, args.out_json)
            sys.exit(0)

        s0_anchor_key = s3_key_exp("S0_no_removal", imp, method, False)
        if s3_exists(s3, s0_anchor_key):
            # Case B: S0 harmonized present — derive missing outputs without Octave
            _log(
                label,
                f"S0 harmonized found on S3 ({len(missing_specs)} output(s) missing). "
                "Deriving without re-harmonizing.",
            )
            try:
                rows = _run_reuse_path(
                    s3=s3,
                    imp=imp,
                    method=method,
                    missing_specs=missing_specs,
                    all_output_specs=all_output_specs,
                    args=args,
                    label=label,
                )
            except Exception:
                _log(
                    label,
                    f"Reuse path failed — falling through to full harmonization:\n"
                    f"{traceback.format_exc()}",
                )
            else:
                _write_sidecar(rows, s3, imp, method, args.out_json)
                sys.exit(0)

        # Case C: S0 harmonized not on S3 — fall through to full harmonization
        _log(label, "S0 harmonized not found on S3 — running full harmonization.")

    # ── Step 1: Memory pre-flight ──────────────────────────────────────────────
    ram = free_gb()
    if ram < args.memory_limit_gb:
        print(
            f"[{_ts()}][{label}] Insufficient RAM: {ram:.1f} GB free < "
            f"{args.memory_limit_gb:.1f} GB required.",
            file=sys.stderr,
        )
        sys.exit(2)
    _log(label, f"RAM pre-flight OK: {ram:.1f} GB free.")

    t_total = time.time()

    # ── Step 2: Download S0 expression and annotation ──────────────────────────
    key_exp = s3_key_prepared_exp("S0_no_removal", imp)
    key_ann = s3_key_prepared_ann("S0_no_removal", imp)

    if not s3_exists(s3, key_exp) or not s3_exists(s3, key_ann):
        print(
            f"[{_ts()}][{label}] S0 prepared dataset not found on S3 "
            f"(key: {key_exp}). Run run_prep_parallel.py first.",
            file=sys.stderr,
        )
        sys.exit(3)

    _log(label, "Downloading S0 expression and annotation ...")
    t0 = time.time()
    exp_full = download_exp_from_s3(s3, key_exp)
    ann_full = download_ann_from_s3(s3, key_ann)
    _log(label, f"S0 ready: {exp_full.shape[0]} samples × {exp_full.shape[1]} genes ({time.time()-t0:.0f}s)")

    # ── Step 3: Download calibration P and Q ───────────────────────────────────
    p_key, q_key = SHAMBHALA_VARIANTS[method]
    _log(label, f"Downloading P={p_key} ({P_DATASETS[p_key]}) and Q={q_key} ({Q_DATASETS[q_key]}) ...")
    t0 = time.time()
    p_df = download_calib_df(s3, P_DATASETS[p_key])
    q_df = download_calib_df(s3, Q_DATASETS[q_key])
    _log(label, f"P: {p_df.shape}, Q: {q_df.shape} ({time.time()-t0:.0f}s)")

    if p_df.isna().any().any():
        n_na_cols = int(p_df.isna().any(axis=0).sum())
        raise ValueError(
            f"P calibration '{P_DATASETS[p_key]}' has NaN in {n_na_cols} gene columns. "
            "Fix the calibration CSV before running the benchmark."
        )
    if q_df.isna().any().any():
        n_na_cols = int(q_df.isna().any(axis=0).sum())
        raise ValueError(
            f"Q calibration '{Q_DATASETS[q_key]}' has NaN in {n_na_cols} gene columns."
        )

    # ── Step 4: Run Shambhala harmonization on full S0 dataset ────────────────
    _log(label, "Starting Shambhala harmonization on full S0 dataset ...")
    t0 = time.time()

    # Gene intersection across input, P, Q
    input_T = exp_full.T
    p_T = p_df.T
    q_T = q_df.T
    common_genes = input_T.index.intersection(p_T.index).intersection(q_T.index)
    n_dropped_input = len(input_T.index) - len(input_T.index.intersection(common_genes))
    n_dropped_p = len(p_T.index) - len(p_T.index.intersection(common_genes))
    n_dropped_q = len(q_T.index) - len(q_T.index.intersection(common_genes))
    coverage_pct = 100 * len(common_genes) / max(len(input_T.index), 1)
    _log(
        label,
        f"Gene intersection: {len(common_genes)} common ({coverage_pct:.1f}% of input) "
        f"(dropped {n_dropped_input} from input, {n_dropped_p} from P, {n_dropped_q} from Q)",
    )
    if len(common_genes) == 0:
        msg = (
            f"Gene intersection is empty for method={method!r}. "
            f"P='{P_DATASETS[p_key]}' ({len(p_T.index)} genes) and/or "
            f"Q='{Q_DATASETS[q_key]}' ({len(q_T.index)} genes) have no gene names in "
            f"common with S0 input ({len(input_T.index)} genes). "
            "Likely cause: calibration CSV stored in genes×samples orientation instead of FL convention."
        )
        print(f"[{_ts()}][{label}] ABORT: {msg}", file=sys.stderr)
        rows = _make_failed_rows(imp, method, status="empty_intersection")
        _write_sidecar(rows, s3, imp, method, args.out_json)
        sys.exit(1)
    input_T = input_T.loc[common_genes]
    p_T = p_T.loc[common_genes]
    q_T = q_T.loc[common_genes]

    # Approach C: centroid-based P subsampling
    if args.max_p_samples is not None and p_T.shape[1] > args.max_p_samples:
        original_p_size = p_T.shape[1]
        p_df_aligned = p_T.T
        p_mean = p_df_aligned.mean(axis=0)
        distances = ((p_df_aligned - p_mean) ** 2).sum(axis=1)
        p_T = p_T[distances.nsmallest(args.max_p_samples).index]
        _log(label, f"P subsampled from {original_p_size} to {args.max_p_samples} samples (centroid).")

    # NA strategy (always 'drop' for pre-imputed S0 data)
    clean_df, na_mask = na_handling.apply_na_strategy(
        input_T.T, strategy=args.na_strategy
    )
    clean_T = clean_df.T
    p_T = p_T.loc[clean_T.index]
    _log(label, f"NA strategy '{args.na_strategy}': {len(na_mask.dropped_genes)} genes dropped.")

    # Approaches A + D: precompute QN reference from P; Python-side QN
    p_for_octave = p_T
    if args.precompute_qn_reference or args.synthetic_cublock_p:
        p_sorted = np.sort(p_T.values, axis=0)
        p_qn_reference = p_sorted.mean(axis=1)

    if args.precompute_qn_reference:
        import qnorm as _qnorm  # noqa: PLC0415

        clean_T = _qnorm.quantile_normalize(clean_T, axis=1, target=p_qn_reference)
        _log(label, f"Python-side QN applied against P reference ({p_T.shape[1]} samples).")

    if args.synthetic_cublock_p:
        p_for_octave = pd.DataFrame(
            p_qn_reference[:, np.newaxis],
            index=p_T.index,
            columns=["P_qn_reference"],
        )
        _log(label, "Synthetic 1-column P will be passed to Octave CuBlock.")

    # Approach E: precompute k-means gene clusters from P
    fixed_clusters = None
    if args.precompute_cublock_clusters:
        from shambhala.octave_bridge import precompute_clusters_from_p  # noqa: PLC0415

        _log(label, f"Precomputing k-means gene clusters from P ({p_for_octave.shape[1]} samples) ...")
        fixed_clusters = precompute_clusters_from_p(
            p_T=p_for_octave,
            k=5,
            random_seed=args.random_seed,
        )
        _log(label, "Gene cluster assignments precomputed.")

    # Q statistics
    q_df_filtered = q_T.loc[clean_T.index].T
    rm, rs = q_rescale.compute_q_statistics(q_df_filtered)
    _log(label, f"Q statistics computed for {len(rm)} genes.")

    try:
        exp_harmonized_T = parallel_module.harmonize_parallel(
            input_df=clean_T,
            p_df=p_for_octave,
            rm=rm,
            rs=rs,
            k=5,
            n_workers=args.n_shambhala_workers,
            octave_scripts_dir=str(_OCTAVE_DIR),
            octave_bin=args.octave_bin,
            timeout_s=args.shambhala_timeout_s,
            random_seed=args.random_seed,
            fixed_clusters=fixed_clusters,
            python_cublock=args.python_cublock,
        )
        # Restore NAs and transpose to FL convention (samples × genes)
        exp_harmonized = na_handling.restore_na_genes(exp_harmonized_T.T, na_mask)
        _log(label, f"Harmonization done: {exp_harmonized.shape} ({time.time()-t0:.0f}s)")
    except Exception:
        print(
            f"[{_ts()}][{label}] Shambhala harmonization FAILED:\n{traceback.format_exc()}",
            file=sys.stderr,
        )
        rows = _make_failed_rows(imp, method, status="failed")
        _write_sidecar(rows, s3, imp, method, args.out_json)
        sys.exit(1)

    # ── Step 5: Free raw S0 expression ────────────────────────────────────────
    del exp_full, ann_full, clean_T, p_T, q_T, q_df_filtered, input_T, clean_df
    gc.collect()
    _log(label, f"Raw S0 freed. RAM: {free_gb():.1f} GB free.")

    # ── Step 6: Download all 14 strategy annotations ───────────────────────────
    _log(label, "Downloading 14 strategy annotations ...")
    t0 = time.time()
    ann_strat: dict[str, pd.DataFrame] = {}
    for strat in ALL_STRATEGIES:
        ann_strat[strat] = download_ann_from_s3(s3, s3_key_prepared_ann(strat, imp))
    counts = {strat: len(ann_strat[strat]) for strat in ALL_STRATEGIES}
    _log(label, f"Annotations ready ({time.time()-t0:.0f}s). Sample counts: {counts}")

    # ── Step 7: Write all 24 (strategy × post_rm) outputs ─────────────────────
    rows: list[dict] = []
    n_done = 0
    for strat in ALL_STRATEGIES:
        for post_rm in [False, True]:
            n_done += 1
            tag = f"{strat}×post_rm={post_rm}"

            # Filter harmonized matrix to strategy sample IDs
            strat_sample_ids = exp_harmonized.index.intersection(ann_strat[strat].index)
            exp_strat = exp_harmonized.loc[strat_sample_ids]
            ann_s = ann_strat[strat].loc[strat_sample_ids]
            _log(label, f"[{n_done}/24] {tag}: {len(ann_s)} samples after filter")

            if post_rm:
                try:
                    outliers = identify_outlier_batches(exp_strat, ann_s)
                    if outliers:
                        _log(label, f"  Removing {len(outliers)} outlier batch(es): {outliers}")
                        ann_s = ann_s[~ann_s[BATCH_COL].isin(outliers)]
                        exp_strat = exp_strat.loc[ann_s.index]
                    else:
                        _log(label, "  No outlier batches found.")
                except Exception:
                    _log(label, f"  Post-removal failed (continuing):\n{traceback.format_exc()}")

            key = s3_key_exp(strat, imp, method, post_rm)

            if args.skip_if_exists and s3_exists(s3, key):
                _log(label, f"  Already on S3 — skipping {tag}.")
                rows.append({
                    "strat": strat, "imp": imp, "method": method,
                    "post_rm": post_rm,
                    "r2_batch": float("nan"), "r2_diag": float("nan"),
                    "n_samples": float("nan"), "n_genes": float("nan"),
                    "status": "cached",
                })
                del exp_strat, ann_s
                gc.collect()
                continue

            try:
                r2_b = r2_batch(exp_strat, ann_s, batch_col=BATCH_COL)
                r2_d = r2_batch(exp_strat, ann_s, batch_col=BIO_COL)
            except Exception:
                r2_b = r2_d = float("nan")

            _log(label, f"  r2_batch={r2_b:.3f} r2_diag={r2_d:.3f} — uploading {tag} ...")
            try:
                upload_exp_to_s3(exp_strat, s3, key)
                status = "ok"
                _log(label, f"  Upload complete.")
            except Exception:
                status = "upload_failed"
                _log(label, f"  Upload FAILED (continuing):\n{traceback.format_exc()}")

            rows.append({
                "strat": strat, "imp": imp, "method": method,
                "post_rm": post_rm,
                "r2_batch": r2_b, "r2_diag": r2_d,
                "n_samples": len(ann_s), "n_genes": exp_strat.shape[1],
                "status": status,
            })
            del exp_strat, ann_s
            gc.collect()

    # ── Step 8: Write and upload sidecar JSON ──────────────────────────────────
    _write_sidecar(rows, s3, imp, method, args.out_json)

    statuses = [r["status"] for r in rows]
    total_elapsed = time.time() - t_total
    _log(
        label,
        f"Done: {statuses.count('ok')} ok, {statuses.count('cached')} cached, "
        f"{statuses.count('upload_failed')} upload_failed, "
        f"{statuses.count('failed')} failed. "
        f"Total elapsed: {total_elapsed:.0f}s",
    )
    sys.exit(0)


if __name__ == "__main__":
    main()
