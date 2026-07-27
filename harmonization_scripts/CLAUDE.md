# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Directory Contains

Cross-product benchmark scripts for the Shambhala2 harmonization method. Runs 3 imputation methods × 18 P/Q calibration variants = **54 worker jobs**, each writing 28 (strategy × post_rm) outputs → **1,512 total S3 objects** (original 12-strategy run produced 1,296; the 216 new objects for J/K are derived via Case B without re-running Octave).

Parent package context: `Shambhala_containerized/CLAUDE.md`.

---

## Commands

All commands must be run from `Shambhala_containerized/` (the package root), not from this subdirectory.

```bash
# Configuration — required. shambhala_bench_shared.py raises KeyError at import without it.
export SHAMBHALA_S3_BUCKET=your-bucket

# Smoke test — single job, strict imputation
python harmonization_scripts/run_shambhala_job.py \
    --imp strict \
    --method shambhala_P0std_Q0std \
    --n-shambhala-workers 2 \
    --octave-bin /usr/bin/octave \
    --random-seed 42 \
    --out-json /tmp/test_result.json

# Full 54-job run (K8s pod; recommended config for 16-vCPU shambhala pod)
nohup python harmonization_scripts/run_shambhala_parallel.py \
    --n-workers 1 --n-shambhala-workers 16 \
    --octave-bin /usr/bin/octave --skip-if-exists \
    --memory-limit-gb 10.0 --timeout-s 7200 --random-seed 42 \
    > /workspace/shambhala_run.log 2>&1 &

# NBGPL570 variant (NP=250, slowest): add --precompute-qn-reference for 2-4× speedup
nohup python harmonization_scripts/run_shambhala_parallel.py \
    --n-workers 1 --n-shambhala-workers 16 \
    --octave-bin /usr/bin/octave --skip-if-exists \
    --precompute-qn-reference \
    --memory-limit-gb 10.0 --timeout-s 14400 --shambhala-timeout-s 7200 --random-seed 42 \
    > /workspace/shambhala_run.log 2>&1 &

# Derive strategy I outputs without re-harmonizing (requires 54 S0 jobs already done)
python harmonization_scripts/derive_rare_batches_outputs.py \
    --skip-if-exists --log-level INFO

# Speed-up test suite (run from Shambhala_containerized/)
pytest tests/speed_up_tests/ -v -s

# Check S3 progress
aws s3 ls s3://your-bucket/FL_batch_correction/exp/ \
    | grep shambhala | wc -l   # target: 1512 (after J/K added); 1296 for original 12-strategy run
```

---

## Architecture

### Three-file pipeline

| File | Role |
|---|---|
| `shambhala_bench_shared.py` | Shared constants, S3 helpers, PCA R² metrics. **Must not import `rpy2`** — copied verbatim from main `bench_shared.py` to avoid its R dependency. |
| `run_shambhala_job.py` | Worker: one `(imp × method)` pair → three-case execution (see below); writes 28 S3 outputs + JSON sidecar. |
| `run_shambhala_parallel.py` | Dispatcher: `ThreadPoolExecutor` over `run_shambhala_job.py` subprocesses; memory guard, failed-jobs log, consolidated `shambhala_metrics.csv`. |
| `derive_rare_batches_outputs.py` | Fast path: derives `I_rare_batches_removed` by filtering rows of existing S0 outputs — no Octave re-run. |

### Why S0 → all strategies (key design insight)

Shambhala normalizes each sample independently. A job runs Shambhala once on the full `S0_no_removal` dataset, then filters post-hoc to each of the 14 strategy sample sets. This collapses 54 × 14 = 756 potential Octave runs to just 54.

### Three-case execution in `run_shambhala_job.py`

Selected automatically when `--skip-if-exists` is set:

| Case | Condition | Action | Octave invoked? |
|------|-----------|--------|-----------------|
| A | All 28 outputs present on S3 | Early exit | No |
| B | `S0_no_removal__post0` present, ≥1 output missing | Download S0 file, filter rows per strategy, upload missing | **No** |
| C | S0 harmonized absent, or `--skip-if-exists` not set | Full harmonization from Octave | **Yes** |

Case B enables adding new batch-removal strategies without re-running Octave: the `S0_no_removal__post0` file already contains all harmonized samples; missing strategy outputs are derived by row-filtering it.

If the S0 download fails in Case B, the job logs a warning and falls through to Case C automatically.

### Exit codes for `run_shambhala_job.py`

| Code | Meaning |
|---|---|
| 0 | Success (including Case A all-cached and Case B reuse) |
| 1 | Harmonization failed |
| 2 | Insufficient RAM at startup |
| 3 | S0 prepared dataset missing from S3 — run `run_prep_parallel.py` first |

### S3 key pattern (identical to main benchmark — compatible with `harmonization-metrics/`)

```
FL_batch_correction/exp/{strat}__{imp}__{method}__post{0|1}.tsv.gz
FL_batch_correction/shambhala_metrics/{imp}__{method}.json
FL_batch_correction/shambhala_metrics.csv
```

### Two-level timeout hierarchy

`--timeout-s` controls how long the dispatcher waits for the entire `run_shambhala_job.py` subprocess. `--shambhala-timeout-s` controls the inner per-Octave-batch limit. Defaults to `--timeout-s` when omitted. Both must be large enough for the batch size `NH = n_S0_samples / n_shambhala_workers`.

With `--n-shambhala-workers 5`, NH ≈ 1435 — this exceeded 7200 s per batch (timed out in the 2026-05-17 run). Use `--n-shambhala-workers 16` (NH ≈ 449) on the 16-vCPU pod.

### Speed-up flags (integrity-ranked)

| Flag | Approach | Max deviation | Notes |
|---|---|---|---|
| `--precompute-qn-reference` | A+D | < 0.05 | Python QN from P; recommended for NBGPL570 (NP=250) |
| `--max-p-samples N` | C | < 0.10 | Centroid subsampling of P |
| `--precompute-cublock-clusters` | E | < 0.10 | sklearn k-means reused per sample |
| `--python-cublock` | G | < 0.10 | Full Python CuBlock; no Octave for CuBlock step |

Always validate a new combination with `pytest tests/speed_up_tests/ -v` before a production run.

### Failed-jobs log

Failed `imp__method` keys are written to `failed_jobs_shambhala.txt` and mirrored to S3 under `FL_batch_correction/failed_shambhala_{POD_NAME}.txt`. Use `--retry-failed` to re-attempt them; use `--skip-if-exists` to skip already-completed ones.

---

## K8s Pod

Pod name: `shambhala-bench`, namespace `<your-namespace>`. Base image `python:3.12-bookworm`, Octave via apt (`/usr/bin/octave`). 16 vCPU / 128 GiB. Full operations guide: `k8s/README.md`.

Sync code from Mac:
```bash
rsync -avz --progress \
    -e "ssh -p 2223 -i ~/.ssh/id_ed25519 -o StrictHostKeyChecking=no" \
    /Users/user890/Desktop/fl_subset/shambhala_adoption/Shambhala_containerized/ \
    root@localhost:/app/Shambhala_containerized/
```

After sync, inside pod: `cd /app/Shambhala_containerized && pip install -e .`

---

## Calibration Dataset Registry

9 P datasets × 2 Q datasets = 18 `SHAMBHALA_VARIANTS` (keys like `shambhala_P0std_Q0std`). Stored on S3 under `FL_batch_correction/calibration_datasets/` in FL convention (rows = samples, columns = genes). All calibration CSVs must be NaN-free — the job aborts with a clear error if NaNs are found.

The slowest variant is `shambhala_NBGPL570_Q0std` (NP=250 vs. NP=39 for `P0std`): QN and CuBlock both scale with NP, causing a ~14–15× slowdown. Always use `--precompute-qn-reference` for this variant.
