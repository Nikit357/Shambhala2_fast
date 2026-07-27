# Shambhala Cross-Product Benchmark

This directory contains scripts to run the Shambhala2 harmonization algorithm across a full
cross-product of imputation methods, calibration datasets, batch-removal strategies, and
post-removal variants.

**Key insight:** Shambhala normalizes each sample independently. Running Shambhala once on
the full S0 dataset covers all 12 batch-removal strategies — the harmonized matrix is filtered
post-hoc to each strategy's sample set.

---

## Configuration

Every script here reads the target S3 bucket from the environment and raises at import time if
it is unset:

```bash
export SHAMBHALA_S3_BUCKET=your-bucket
```

See `.env.example` in the repository root. There is deliberately no default — a wrong-bucket
default would fail silently or write somewhere unintended.

---

## What this produces

| Dimension | Values | Count |
|---|---|---|
| Imputation | `strict`, `knn`, `softimpute` | 3 |
| Shambhala variant (P × Q) | 18 combinations | 18 |
| Batch-removal strategy | 14 strategies (same as main framework) | 14 |
| Post-removal | `False` (post0), `True` (post1) | 2 |
| **Total S3 outputs** | 3 × 18 × 14 × 2 | **1,512** |
| **Total Shambhala runs** | 3 × 18 | **54** |

Output S3 keys follow the **same pattern as the main benchmark** — compatible with
`harmonization-metrics/` and `harmonization_benchmark.ipynb` without any changes.

---

## Prerequisites

- Pod running (see `k8s/README.md`).
- Code synced to pod via rsync.
- `pip install -e .` done inside the pod (installs the `shambhala` package).
- S3 credentials mounted at `/root/.aws/`.
- All 3 S0 prepared datasets present on S3:
  `FL_batch_correction/prepared/S0_no_removal__{strict,knn,softimpute}__exp.tsv.gz`

---

## Quick start

### Full run (all 54 jobs)

```bash
nohup python harmonization_scripts/run_shambhala_parallel.py \
    --n-workers 3 \
    --n-shambhala-workers 5 \
    --octave-bin /usr/bin/octave \
    --skip-if-exists \
    --memory-limit-gb 10.0 \
    --timeout-s 21600 \
    --random-seed 42 \
    > /workspace/shambhala_run.log 2>&1 &
```

### Smoke test (1 variant, strict imputation only)

```bash
python harmonization_scripts/run_shambhala_job.py \
    --imp strict \
    --method shambhala_P0std_Q0std \
    --n-shambhala-workers 2 \
    --octave-bin /usr/bin/octave \
    --random-seed 42 \
    --out-json /tmp/test_result.json
```

### Targeted subset (2 variants, strict only)

```bash
python harmonization_scripts/run_shambhala_parallel.py \
    --imps strict \
    --methods shambhala_P0std_Q0std,shambhala_NBBags_Q0std \
    --n-workers 2 \
    --n-shambhala-workers 3 \
    --skip-if-exists \
    --octave-bin /usr/bin/octave
```
### Last real world example of Daniil's run in May 20th
```bash
python harmonization_scripts/run_shambhala_parallel.py   --imps strict    --methods shambhala_P0std_QNBKass,shambhala_ANTE_QNBKass,shambhala_NBGPL570_QNBKass,shambhala_NBKass_QNBKass,shambhala_NBRNAseq_QNBKass   --n-workers 1  --n-shambhala-workers 30   --skip-if-exists  --octave-bin /usr/bin/octave > shambhala_run_250620.log
```

### Derive `I_rare_batches_removed` from existing S0 outputs (no re-harmonization)

If the 54-job run has already completed, strategy I outputs can be derived cheaply
by filtering rows from the existing S0 harmonized matrices:

```bash
python harmonization_scripts/derive_rare_batches_outputs.py \
    --skip-if-exists \
    --log-level INFO \
    > /workspace/derive_rare_batches.log 2>&1
```

---

## Full command reference

### `run_shambhala_parallel.py` (dispatcher)

| Flag | Default | Description |
|---|---|---|
| `--n-workers` | 3 | Parallel `run_shambhala_job.py` subprocesses. |
| `--n-shambhala-workers` | 5 | Octave parallelism inside each job. |
| `--octave-bin` | `octave` | Path to Octave binary. Always pass `/usr/bin/octave` inside pod. |
| `--memory-limit-gb` | 10.0 | Wait until free RAM ≥ this before launching each job. |
| `--timeout-s` | 21600 | Per-job timeout (6 h). |
| `--skip-if-exists` | — | Skip jobs where all 28 S3 outputs already exist. |
| `--retry-failed` | — | Re-attempt jobs from `failed_jobs_shambhala.txt`. |
| `--imps` | all | Comma-separated imputation subset (e.g. `strict,knn`). |
| `--methods` | all | Comma-separated method subset. |
| `--tmp-dir` | `/tmp/shambhala_jobs` | Directory for per-job JSON sidecars. |
| `--random-seed` | None | Random seed for reproducible k-means in Octave. |

### `run_shambhala_job.py` (worker)

| Flag | Default | Description |
|---|---|---|
| `--imp` | required | Imputation method. |
| `--method` | required | Shambhala variant key. |
| `--skip-if-exists` | — | Skip per-output keys already on S3. |
| `--out-json` | None | Local path for JSON sidecar. |
| `--memory-limit-gb` | 8.0 | Abort if free RAM < this at startup. |
| `--n-shambhala-workers` | 5 | Octave parallel workers (capped at 10). |
| `--octave-bin` | `octave` | Path to Octave binary. |
| `--random-seed` | None | Seed for Octave k-means initialization. |
| `--shambhala-timeout-s` | 7200 | Per-batch Octave subprocess timeout. |
| `--na-strategy` | `drop` | NA handling: `drop` (default) or `knn`. Always use `drop` for pre-imputed S0 data. |

---

## 18-method variant table

| Method key | P dataset | Q dataset |
|---|---|---|
| `shambhala_P0std_Q0std` | P0_standard.csv | Q0_standard.csv |
| `shambhala_ANTE_Q0std` | ANTE.csv | Q0_standard.csv |
| `shambhala_GTExAffy_Q0std` | GTExAffymetrix.csv | Q0_standard.csv |
| `shambhala_NBBags_Q0std` | Normal_B_BG_BAGS.csv | Q0_standard.csv |
| `shambhala_NBGPL570_Q0std` | Normal_B_GPL570.csv | Q0_standard.csv |
| `shambhala_NBKass_Q0std` | Normal_B_Kassandra.csv | Q0_standard.csv |
| `shambhala_NBRNAseq_Q0std` | Normal_B_RNASeq.csv | Q0_standard.csv |
| `shambhala_NBext_Q0std` | Normal_B_cells_extended.csv | Q0_standard.csv |
| `shambhala_Oncobox_Q0std` | OncoboxCancer.csv | Q0_standard.csv |
| `shambhala_P0std_QNBKass` | P0_standard.csv | Normal_B_Kassandra.csv |
| `shambhala_ANTE_QNBKass` | ANTE.csv | Normal_B_Kassandra.csv |
| `shambhala_GTExAffy_QNBKass` | GTExAffymetrix.csv | Normal_B_Kassandra.csv |
| `shambhala_NBBags_QNBKass` | Normal_B_BG_BAGS.csv | Normal_B_Kassandra.csv |
| `shambhala_NBGPL570_QNBKass` | Normal_B_GPL570.csv | Normal_B_Kassandra.csv |
| `shambhala_NBKass_QNBKass` | Normal_B_Kassandra.csv | Normal_B_Kassandra.csv |
| `shambhala_NBRNAseq_QNBKass` | Normal_B_RNASeq.csv | Normal_B_Kassandra.csv |
| `shambhala_NBext_QNBKass` | Normal_B_cells_extended.csv | Normal_B_Kassandra.csv |
| `shambhala_Oncobox_QNBKass` | OncoboxCancer.csv | Normal_B_Kassandra.csv |

**Note on `shambhala_NBKass_QNBKass`:** P and Q are the same file (`Normal_B_Kassandra.csv`).
P anchors quantile normalization shape and Q defines the target mean/std. Using the same file
creates a circular normalization target but is numerically valid.

---

## Timing estimates

| Imputation | Approx. time per variant | Notes |
|---|---|---|
| `strict` | 30–60 min | Small gene set (3,447 genes), fewer NAs |
| `knn` | 2–4 h | Larger gene set, ~7,238 samples |
| `softimpute` | 2–4 h | Similar to knn |

Total wall time with `--n-workers 3`: strict finishes in ~3–4 h; knn/softimpute take 12–20 h.

---

## Resuming a partial run

```bash
# Skip completed jobs (checks S3)
python harmonization_scripts/run_shambhala_parallel.py --n-workers 3 --skip-if-exists ...

# Retry previously failed jobs
python harmonization_scripts/run_shambhala_parallel.py --n-workers 3 --retry-failed --skip-if-exists ...
```

Failed jobs are logged in `harmonization_scripts/failed_jobs_shambhala.txt` (one `imp__method` key per line)
and mirrored to `s3://your-bucket/FL_batch_correction/failed_shambhala_{POD_NAME}.txt`.

---

## Monitoring

```bash
tail -f /workspace/shambhala_run.log

aws s3 ls s3://your-bucket/FL_batch_correction/exp/ \
    | grep shambhala | wc -l    # target: 1512 (after J/K added); 1296 for original 12-strategy run

htop    # CPU / RAM usage
```

---

## Output format

All outputs are gzip-compressed TSV files:
- Rows = samples (index column = sample ID).
- Columns = HGNC gene symbols.
- Values = log-transformed harmonized expression.

This is the same format as the main benchmark. Load with:

```python
import pandas as pd
df = pd.read_csv("s3://...", sep="\t", index_col=0, compression="gzip")
```

---

## NA strategy note

Always use `--na-strategy drop` (the default). The S0 `knn` and `softimpute` prepared datasets
are already imputed — double-imputing would distort values. The `drop` strategy removes any
residual NaN genes before Octave and marks them as NaN in the output.

---

## S3 output paths

```
s3://your-bucket/FL_batch_correction/
    exp/{strat}__{imp}__{method}__post{0|1}.tsv.gz    # harmonized expression (×1,512)
    shambhala_metrics/{imp}__{method}.json             # per-job sidecar (×54)
    shambhala_metrics.csv                              # consolidated metrics (all rows)
```

---

## Performance Tuning

### Why NBGPL570 is slow

Both QN and CuBlock inside each Octave subprocess scale with the number of P samples (NP).
`shambhala_P0std_Q0std` (NP=39) runs at ~10 s/sample; `shambhala_NBGPL570_Q0std` (NP=250)
takes ~155 s/sample — a 14–15× slowdown for only 6.4× more P samples.

### Flag priority ranking (highest integrity first)

| Priority | Flag | Approach | Max deviation | Notes |
|---|---|---|---|---|
| 1 (always active) | *(built-in)* | **B** | 0 | Skip wasted P-column CuBlock fits |
| 2 | `--precompute-qn-reference` | **A+D** | < 0.05 | Python QN from P; recommended for NBGPL570 |
| 3 | `--max-p-samples N` | **C** | < 0.10 | Centroid P subsampling |
| 4 | `--precompute-cublock-clusters` | **E** | < 0.10 | sklearn k-means clusters reused per sample |
| 5 | `--python-cublock` | **G** | < 0.10 | Full Python CuBlock; no Octave for CuBlock |

### Three worked example commands (K8s pod)

**1. Safe/fast — recommended for production NBGPL570 runs:**
```bash
nohup python harmonization_scripts/run_shambhala_parallel.py \
    --n-workers 3 --n-shambhala-workers 10 \
    --octave-bin /usr/bin/octave --skip-if-exists \
    --precompute-qn-reference \
    --memory-limit-gb 10.0 --timeout-s 21600 --random-seed 42 \
    > /workspace/shambhala_run.log 2>&1 &
```

**2. Maximum speed — use when time matters more than exact reproducibility:**
```bash
nohup python harmonization_scripts/run_shambhala_parallel.py \
    --n-workers 3 --n-shambhala-workers 15 \
    --octave-bin /usr/bin/octave --skip-if-exists \
    --precompute-qn-reference --synthetic-cublock-p \
    --max-p-samples 40 --precompute-cublock-clusters \
    --memory-limit-gb 10.0 --timeout-s 21600 --random-seed 42 \
    > /workspace/shambhala_fast.log 2>&1 &
```

**3. Full Python — no Octave for CuBlock (experimental validation):**
```bash
python harmonization_scripts/run_shambhala_parallel.py \
    --n-workers 1 --n-shambhala-workers 1 \
    --octave-bin /usr/bin/octave \
    --precompute-qn-reference --python-cublock \
    --random-seed 42 \
    --imps strict --methods shambhala_P0std_Q0std
```

### Validation

Before switching to a faster combination in production, run the integration tests
against the expected_output fixture:

```bash
cd /app/Shambhala_containerized
pytest tests/speed_up_tests/ -v -s 2>&1 | tee /workspace/speed_up_test_results.txt
```

Check that `max_abs_diff` for every combination is within the declared tolerance
(see `tests/speed_up_tests/test_combinations.py` for the full table).
