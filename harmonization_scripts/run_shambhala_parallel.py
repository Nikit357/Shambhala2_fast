"""
Shambhala cross-product dispatcher: enumerate all 54 (imp × method) jobs, launch
run_shambhala_job.py subprocesses, collect results, upload consolidated metrics.

3 imputations × 18 Shambhala variants = 54 worker jobs
Each worker writes 28 (strategy × post_rm) outputs → 1,512 total S3 objects.

Usage
-----
# Full run (all 54 jobs):
nohup python harmonization_scripts/run_shambhala_parallel.py \\
    --n-workers 3 \\
    --n-shambhala-workers 5 \\
    --octave-bin /usr/bin/octave \\
    --skip-if-exists \\
    --memory-limit-gb 10.0 \\
    --timeout-s 21600 \\
    --random-seed 42 \\
    > /workspace/shambhala_run.log 2>&1 &

# Targeted subset:
python harmonization_scripts/run_shambhala_parallel.py \\
    --imps strict \\
    --methods shambhala_P0std_Q0std,shambhala_NBBags_Q0std \\
    --n-workers 2 --n-shambhala-workers 3 --skip-if-exists \\
    --octave-bin /usr/bin/octave

# Resume after failure:
python harmonization_scripts/run_shambhala_parallel.py \\
    --n-workers 3 --skip-if-exists --retry-failed

--n-workers N
    Number of parallel run_shambhala_job.py subprocesses (default: 3).
    Each job runs up to --n-shambhala-workers Octave processes internally.

--n-shambhala-workers N
    Internal Octave parallelism passed through to each job (default: 5).
    Total Octave processes ≤ n-workers × n-shambhala-workers.

--memory-limit-gb F
    Dispatcher pauses until free RAM ≥ F GB before launching each job (default: 10.0).
    Each worker receives F/2 as its own --memory-limit-gb.

--skip-if-exists
    At startup, scans S3 for existing outputs. Jobs with all 28 outputs already
    present are marked cached and not relaunched.

--retry-failed
    Re-attempt jobs listed in failed_jobs_shambhala.txt.

--timeout-s N
    Per-job timeout in seconds (default: 21600 = 6 hours). Controls how long
    the dispatcher waits for the entire run_shambhala_job.py subprocess.

--shambhala-timeout-s N
    Per-Octave-batch timeout in seconds passed through to each job subprocess as
    --shambhala-timeout-s. This is the limit applied to each individual Octave
    call inside the job. Defaults to --timeout-s if not set.
    With --n-shambhala-workers 16 (NH≈449 per batch), 7200 s is a safe value.
    With --n-shambhala-workers 5  (NH≈1435 per batch), even 7200 s is not
    enough — Octave times out. Increase workers to reduce NH before lowering
    this timeout.

--imps LIST
    Comma-separated imputation subset (default: strict,knn,softimpute).

--methods LIST
    Comma-separated method subset (default: all 18).

--tmp-dir PATH
    Directory for per-job JSON sidecars (default: /tmp/shambhala_jobs).

Speed-Up Flags
--------------
All flags below are passed through verbatim to each run_shambhala_job.py subprocess.
See README.md "Speed-Up Options" for integrity ratings and expected speedups.

--precompute-qn-reference
    Approach A+D: compute QN reference from P once in Python; apply via qnorm per sample.
    Minor numerical deviation (max_abs_diff < 0.05). Recommended for NBGPL570.

--synthetic-cublock-p
    Approach A extra: pass a single synthetic P column to CuBlock (mean of QN reference).
    Only useful combined with --precompute-qn-reference.

--max-p-samples N
    Approach C: subsample P to N samples closest to P centroid before QN/CuBlock.

--precompute-cublock-clusters
    Approach E: precompute k-means gene clusters from P via sklearn; inject into
    CuBlock_fixed.m. Eliminates per-sample k-means inside Octave.

--python-cublock
    Approach G: replace Octave CuBlock with pure Python (sklearn + numpy).
    Fastest option; bypasses Octave for CuBlock. Requires qnorm for QN step.
"""
from __future__ import annotations

import argparse
import json
import os
import signal
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import boto3
import botocore.exceptions
import pandas as pd

# Allow running from the harmonization_scripts/ directory
_THIS_DIR = Path(__file__).resolve().parent
_SHAMBHALA_ROOT = _THIS_DIR.parent
if str(_SHAMBHALA_ROOT) not in sys.path:
    sys.path.insert(0, str(_SHAMBHALA_ROOT))

from harmonization_scripts.shambhala_bench_shared import (
    ALL_IMPUTATION,
    ALL_METHODS,
    ALL_STRATEGIES,
    S3_BUCKET,
    S3_PREFIX,
    free_gb,
    s3_key_exp,
    s3_key_sidecar,
    s3_key_shambhala_metrics_csv,
)

_MEM_POLL_S = 30.0
FAILED_LOG_PATH: Path = _THIS_DIR / "failed_jobs_shambhala.txt"
_POD_NAME: str = os.environ.get("POD_NAME", "local")
_S3_FAILED_PREFIX = f"{S3_PREFIX}/failed_shambhala_"


def _ts() -> str:
    return time.strftime("%H:%M:%S")


def _job_key(imp: str, method: str) -> str:
    return f"{imp}__{method}"


# ── Failed-jobs log ────────────────────────────────────────────────────────────

def _load_failed_log() -> set[str]:
    if not FAILED_LOG_PATH.exists():
        return set()
    with open(FAILED_LOG_PATH) as fh:
        return {line.strip() for line in fh if line.strip()}


def _save_failed_log(keys: set[str]) -> None:
    with open(FAILED_LOG_PATH, "w") as fh:
        for key in sorted(keys):
            fh.write(key + "\n")
    print(
        f"Shambhala failed-jobs log written: {FAILED_LOG_PATH} ({len(keys)} entries)",
        flush=True,
    )


def _sync_failed_log_from_s3() -> None:
    s3 = boto3.client("s3")
    try:
        paginator = s3.get_paginator("list_objects_v2")
        s3_keys: set[str] = set()
        for page in paginator.paginate(Bucket=S3_BUCKET, Prefix=_S3_FAILED_PREFIX):
            for obj in page.get("Contents", []):
                body = (
                    s3.get_object(Bucket=S3_BUCKET, Key=obj["Key"])["Body"]
                    .read()
                    .decode()
                )
                s3_keys |= {ln.strip() for ln in body.splitlines() if ln.strip()}
        if s3_keys:
            local_keys = _load_failed_log()
            merged = s3_keys | local_keys
            if merged != local_keys:
                _save_failed_log(merged)
                print(
                    f"[S3 sync] Merged {len(s3_keys)} Shambhala failed keys → {len(merged)} total.",
                    flush=True,
                )
    except botocore.exceptions.ClientError as exc:
        print(f"[S3 sync] Could not read Shambhala failed_jobs from S3: {exc}", flush=True)


def _sync_failed_log_to_s3(keys: set[str]) -> None:
    s3 = boto3.client("s3")
    key = f"{_S3_FAILED_PREFIX}{_POD_NAME}.txt"
    if keys:
        body = "\n".join(sorted(keys)).encode()
        s3.put_object(Bucket=S3_BUCKET, Key=key, Body=body)
        print(
            f"[S3 sync] Uploaded {len(keys)} Shambhala failed keys to s3://{S3_BUCKET}/{key}",
            flush=True,
        )
    else:
        try:
            s3.delete_object(Bucket=S3_BUCKET, Key=key)
        except Exception:
            pass


# ── S3 existence scan ──────────────────────────────────────────────────────────

def _scan_existing_outputs() -> set[str]:
    """
    List all existing exp/ S3 objects matching the shambhala_ pattern.
    Returns a set of full S3 keys present on S3.
    """
    s3 = boto3.client("s3")
    prefix = f"{S3_PREFIX}/exp/"
    existing: set[str] = set()
    paginator = s3.get_paginator("list_objects_v2")
    for page in paginator.paginate(Bucket=S3_BUCKET, Prefix=prefix):
        for obj in page.get("Contents", []):
            key = obj["Key"]
            if "shambhala_" in key:
                existing.add(key)
    return existing


def _all_outputs_exist(imp: str, method: str, existing: set[str]) -> bool:
    """Return True if all 28 output keys for this job already exist on S3."""
    return all(
        s3_key_exp(strat, imp, method, post_rm) in existing
        for strat in ALL_STRATEGIES
        for post_rm in [False, True]
    )


# ── Memory guard ───────────────────────────────────────────────────────────────

def _wait_for_memory(min_free_gb: float) -> None:
    while True:
        if free_gb() >= min_free_gb:
            return
        print(
            f"  [mem-guard] {free_gb():.1f} GB free — pausing {_MEM_POLL_S:.0f}s",
            flush=True,
        )
        time.sleep(_MEM_POLL_S)


# ── Stream forwarding ──────────────────────────────────────────────────────────

def _stream_and_forward(pipe: object, prefix: str) -> None:
    for line in pipe:  # type: ignore[attr-defined]
        print(f"  [{prefix}] {line}", end="", flush=True)


# ── Dispatcher ─────────────────────────────────────────────────────────────────

def _run_dispatcher(
    jobs: list[tuple[str, str]],
    n_workers: int,
    skip_if_exists: bool,
    tmp_dir: Path,
    skip_keys: set[str],
    timeout_s: int,
    shambhala_timeout_s: int,
    memory_limit_gb: float,
    n_shambhala_workers: int,
    octave_bin: str,
    random_seed: int | None,
    precompute_qn_reference: bool = False,
    synthetic_cublock_p: bool = False,
    max_p_samples: int | None = None,
    precompute_cublock_clusters: bool = False,
    python_cublock: bool = False,
) -> tuple[list[dict], set[str], set[str]]:
    jobs_to_run = [(imp, method) for imp, method in jobs if _job_key(imp, method) not in skip_keys]
    skipped_count = len(jobs) - len(jobs_to_run)
    if skipped_count:
        print(f"Skipping {skipped_count} previously failed Shambhala job(s).", flush=True)

    total = len(jobs_to_run)
    print(
        f"Shambhala jobs to run: {total}  |  Workers: {n_workers}  |  "
        f"Shambhala-workers/job: {n_shambhala_workers}  |  "
        f"Timeout (job): {timeout_s}s  |  Timeout (Octave batch): {shambhala_timeout_s}s  |  "
        f"Memory guard: {memory_limit_gb:.1f} GB",
        flush=True,
    )

    def launch(imp: str, method: str) -> tuple[list[dict], float, int]:
        _wait_for_memory(memory_limit_gb)
        tag = f"{imp}×{method}"
        t0 = time.time()
        print(f"[{_ts()}] START {tag} | free={free_gb():.1f}GB", flush=True)

        json_path = tmp_dir / f"{imp}__{method}.json"
        cmd = [
            sys.executable,
            str(_THIS_DIR / "run_shambhala_job.py"),
            "--imp", imp,
            "--method", method,
            "--out-json", str(json_path),
            "--memory-limit-gb", str(memory_limit_gb / 2),
            "--n-shambhala-workers", str(n_shambhala_workers),
            "--octave-bin", octave_bin,
        ]
        if skip_if_exists:
            cmd.append("--skip-if-exists")
        if random_seed is not None:
            cmd.extend(["--random-seed", str(random_seed)])
        cmd.extend(["--shambhala-timeout-s", str(shambhala_timeout_s)])
        if precompute_qn_reference:
            cmd.append("--precompute-qn-reference")
        if synthetic_cublock_p:
            cmd.append("--synthetic-cublock-p")
        if max_p_samples is not None:
            cmd.extend(["--max-p-samples", str(max_p_samples)])
        if precompute_cublock_clusters:
            cmd.append("--precompute-cublock-clusters")
        if python_cublock:
            cmd.append("--python-cublock")

        try:
            proc = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                cwd=str(_SHAMBHALA_ROOT),
            )
            reader = threading.Thread(
                target=_stream_and_forward,
                args=(proc.stdout, tag),
                daemon=True,
            )
            reader.start()
            try:
                proc.wait(timeout=timeout_s)
            except subprocess.TimeoutExpired:
                proc.kill()
                reader.join(timeout=5)
                elapsed = time.time() - t0
                print(f"[{_ts()}] TIMEOUT {tag} after {elapsed:.0f}s", flush=True)
                return ([], elapsed, -9)
            reader.join(timeout=5)
        except Exception as exc:
            elapsed = time.time() - t0
            print(f"[{_ts()}] ERROR {tag}: {exc}", flush=True)
            return ([], elapsed, -1)

        elapsed = time.time() - t0
        rows: list[dict] = []
        if proc.returncode == 0 and json_path.exists():
            try:
                with open(json_path) as fh:
                    rows = json.load(fh)
            except Exception:
                pass
        elif proc.returncode == 3:
            print(f"[{_ts()}] MISSING S0 prepared dataset for {tag}", flush=True)
        elif proc.returncode != 0:
            print(f"[{_ts()}] FAILED {tag} exit={proc.returncode} ({elapsed:.0f}s)", flush=True)
        return (rows, elapsed, proc.returncode)

    all_rows: list[dict] = []
    newly_failed: set[str] = set()
    newly_succeeded: set[str] = set()
    done = 0

    with ThreadPoolExecutor(max_workers=n_workers) as pool:
        future_to_job = {
            pool.submit(launch, imp, method): (imp, method)
            for imp, method in jobs_to_run
        }
        for future in as_completed(future_to_job):
            imp, method = future_to_job[future]
            key = _job_key(imp, method)
            done += 1
            rows, elapsed, rc = future.result()
            if rc == 3:
                # Missing prepared dataset — not a failure, just skip
                continue
            statuses = {r.get("status", "failed") for r in rows} if rows else set()
            if rc != 0 or not rows or statuses & {"failed", "upload_failed"}:
                newly_failed.add(key)
                print(f"  [{done}/{total}] FAILED {imp}×{method} ({elapsed:.0f}s)", flush=True)
            elif statuses <= {"cached"}:
                newly_succeeded.add(key)
                all_rows.extend(rows)
                print(f"  [{done}/{total}] CACHED {imp}×{method} ({elapsed:.0f}s)", flush=True)
            else:
                newly_succeeded.add(key)
                all_rows.extend(rows)
                print(
                    f"  [{done}/{total}] {imp}×{method} ({elapsed:.0f}s) "
                    f"→ {[r['status'] for r in rows[:3]]} ...",
                    flush=True,
                )

    return all_rows, newly_failed, newly_succeeded


def _upload_consolidated_metrics(all_rows: list[dict], s3_sidecar_keys: list[str]) -> None:
    """Collect all per-job JSON sidecars + in-memory rows → shambhala_metrics.csv."""
    s3 = boto3.client("s3")
    import io

    # Try to read any sidecars that weren't in all_rows (from previous partial runs)
    extra_rows: list[dict] = []
    for key in s3_sidecar_keys:
        try:
            body = s3.get_object(Bucket=S3_BUCKET, Key=key)["Body"].read()
            extra_rows.extend(json.loads(body))
        except Exception:
            pass

    combined = all_rows + extra_rows
    if not combined:
        print("No metric rows collected — shambhala_metrics.csv not uploaded.", flush=True)
        return

    df = pd.DataFrame(combined).drop_duplicates(
        subset=["strat", "imp", "method", "post_rm"], keep="last"
    )
    csv_bytes = df.to_csv(index=False).encode()
    buf = io.BytesIO(csv_bytes)
    key = s3_key_shambhala_metrics_csv()
    s3.upload_fileobj(buf, S3_BUCKET, key)
    print(f"Consolidated metrics uploaded: s3://{S3_BUCKET}/{key} ({len(df)} rows)", flush=True)


def main() -> None:
    signal.signal(
        signal.SIGTERM,
        lambda signum, frame: (
            print("\n[dispatcher] SIGTERM — syncing failed-jobs log to S3.", flush=True),
            _sync_failed_log_to_s3(_load_failed_log()),
            sys.exit(0),
        )[-1],
    )

    parser = argparse.ArgumentParser(
        description="Shambhala cross-product dispatcher: 54 jobs → 1,512 S3 outputs."
    )
    parser.add_argument("--n-workers", type=int, default=3)
    parser.add_argument("--skip-if-exists", action="store_true", dest="skip_if_exists")
    parser.add_argument("--retry-failed", action="store_true", dest="retry_failed")
    parser.add_argument("--imps", default=None)
    parser.add_argument("--methods", default=None)
    parser.add_argument("--n-shambhala-workers", type=int, default=5, dest="n_shambhala_workers")
    parser.add_argument("--octave-bin", default="octave", dest="octave_bin")
    parser.add_argument("--memory-limit-gb", type=float, default=10.0, dest="memory_limit_gb")
    parser.add_argument("--timeout-s", type=int, default=21600, dest="timeout_s")
    parser.add_argument(
        "--shambhala-timeout-s",
        type=int,
        default=None,
        dest="shambhala_timeout_s",
        help=(
            "Per-Octave-batch timeout in seconds passed to run_shambhala_job.py "
            "--shambhala-timeout-s. Defaults to --timeout-s if not set. "
            "With --n-shambhala-workers 16 (NH≈449), 7200 s is sufficient. "
            "With --n-shambhala-workers 5 (NH≈1435), 7200 s is not enough."
        ),
    )
    parser.add_argument("--tmp-dir", default="/tmp/shambhala_jobs", dest="tmp_dir")
    parser.add_argument("--random-seed", type=int, default=None, dest="random_seed")
    parser.add_argument(
        "--precompute-qn-reference", action="store_true", default=False,
        dest="precompute_qn_reference",
        help="Precompute QN reference from P; apply Python-side QN per sample (~10-20× speedup).",
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
    args = parser.parse_args()

    imps = args.imps.split(",") if args.imps else ALL_IMPUTATION
    methods = args.methods.split(",") if args.methods else ALL_METHODS
    tmp_dir = Path(args.tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    _sync_failed_log_from_s3()
    previously_failed = _load_failed_log()
    skip_keys = set() if args.retry_failed else previously_failed

    # S3 existence scan
    existing_outputs: set[str] = set()
    if args.skip_if_exists:
        print("[dispatcher] Scanning S3 for existing Shambhala outputs ...", flush=True)
        existing_outputs = _scan_existing_outputs()
        print(f"[dispatcher] Found {len(existing_outputs)} existing Shambhala output files.", flush=True)

    # Build job list
    all_jobs = [(imp, method) for imp in imps for method in methods]
    if args.skip_if_exists:
        cached = [(imp, m) for imp, m in all_jobs if _all_outputs_exist(imp, m, existing_outputs)]
        if cached:
            print(f"[dispatcher] {len(cached)} job(s) fully cached on S3 — will skip.", flush=True)
        # Still pass to dispatcher; worker will do per-output checks for partially-cached jobs
        jobs = all_jobs
    else:
        jobs = all_jobs

    print(
        f"[dispatcher] Total jobs: {len(jobs)} "
        f"({len(imps)} imps × {len(methods)} methods)",
        flush=True,
    )

    all_rows, newly_failed, newly_succeeded = _run_dispatcher(
        jobs=jobs,
        n_workers=args.n_workers,
        skip_if_exists=args.skip_if_exists,
        tmp_dir=tmp_dir,
        skip_keys=skip_keys,
        timeout_s=args.timeout_s,
        shambhala_timeout_s=args.shambhala_timeout_s or args.timeout_s,
        memory_limit_gb=args.memory_limit_gb,
        n_shambhala_workers=args.n_shambhala_workers,
        octave_bin=args.octave_bin,
        random_seed=args.random_seed,
        precompute_qn_reference=args.precompute_qn_reference,
        synthetic_cublock_p=args.synthetic_cublock_p,
        max_p_samples=args.max_p_samples,
        precompute_cublock_clusters=args.precompute_cublock_clusters,
        python_cublock=args.python_cublock,
    )

    updated_failed = (previously_failed - newly_succeeded) | newly_failed
    if updated_failed != previously_failed or newly_failed:
        _save_failed_log(updated_failed)
    _sync_failed_log_to_s3(updated_failed)

    # Collect all sidecar keys to build consolidated metrics
    sidecar_keys = [
        s3_key_sidecar(imp, method)
        for imp in imps
        for method in methods
    ]
    _upload_consolidated_metrics(all_rows, sidecar_keys)

    n_failed = len(newly_failed)
    n_ok = len(newly_succeeded)
    print(
        f"\n[dispatcher] Summary: {n_ok} succeeded, {n_failed} failed. "
        f"Failed keys saved to {FAILED_LOG_PATH}",
        flush=True,
    )


if __name__ == "__main__":
    main()
