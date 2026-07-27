# Cross-Product Shambhala Harmonization — Implementation Plan

**Date:** 2026-05-15  
**Author:** Daniil Nikitin  
**Status:** Approved — ready for implementation

---

## 1. Goal and Motivation

The `Shambhala_containerized/` pipeline currently supports single-run harmonization via `run_shambhala.py`. This plan extends it into a **full cross-product benchmark** matching the structure of the main `harmonization-scripts/` framework but exploiting a key algorithmic property of Shambhala2:

> **Shambhala normalizes each sample independently.** The result for sample *i* does not depend on which other samples are in the input matrix.

This means the 11 sample-removal batch strategies do **not** require 11 separate Shambhala runs. Instead:
1. Run Shambhala once on the full S0 (no removal) dataset per `(imputation × P/Q variant)` combination.
2. Filter the resulting harmonized matrix to each strategy's sample set post-hoc.
3. Apply post-removal (optional outlier batch removal) to each filtered matrix.
4. Upload all results to S3.

The result is a **1,188-output S3 surface** from only **54 Shambhala runs**.

---

## 2. Cross-Product Dimensions

| Dimension | Values | Count |
|---|---|---|
| **Imputation** (S0 prepared datasets) | `strict`, `knn`, `softimpute` | 3 |
| **Shambhala variant** (P × Q calibration) | 18 combinations (see §3) | 18 |
| **Batch-removal strategy** (post-hoc sample filtering) | 11 strategies (same as main framework) | 11 |
| **Post-removal** (outlier batch removal) | `False` (post0), `True` (post1) | 2 |
| **Total S3 outputs** | 3 × 18 × 11 × 2 | **1,188** |
| **Total Shambhala runs** (worker jobs) | 3 × 18 | **54** |
| **Outputs per worker job** | 11 × 2 | **22** |

**Note:** `missforest` is excluded from imputation because the S0 missforest prepared dataset was never successfully produced (>24h/240 GB, confirmed failed in main benchmark). Only `strict`, `knn`, and `softimpute` are available for S0.

---

## 3. Shambhala Variant Registry (18 Methods)

### 3.1 Calibration Datasets

All files are in `s3://your-bucket/FL_batch_correction/calibration_datasets/`.

**Confirmed properties (verified 2026-05-15):**
- All files use **FL convention**: rows = samples, columns = genes. No transposition needed on load.
- All gene columns use **HGNC symbols**. No probe-to-gene mapping step needed.

**Q datasets (2):**

| File | Short key | Description |
|---|---|---|
| `Q0_standard.csv` | `Q0std` | Original Shambhala2 Q reference (mixed tissue RNA-seq from GTEx) |
| `Normal_B_Kassandra.csv` | `QNBKass` | Normal B-cell reference from BG Kassandra project |

**P datasets (9):**

| File | Short key | Description |
|---|---|---|
| `P0_standard.csv` | `P0std` | Original Shambhala2 P reference (39 Affymetrix GPL570 healthy tissue samples) |
| `ANTE.csv` | `ANTE` | ANTE calibration set |
| `GTExAffymetrix.csv` | `GTExAffy` | GTEx samples profiled on Affymetrix |
| `Normal_B_BG_BAGS.csv` | `NBBags` | Normal B cells from BG-BAGS classifier dataset |
| `Normal_B_GPL570.csv` | `NBGPL570` | Normal B cells on GPL570 platform |
| `Normal_B_Kassandra.csv` | `NBKass` | Normal B cells from BG Kassandra project |
| `Normal_B_RNASeq.csv` | `NBRNAseq` | Normal B cells from RNA-seq cohorts |
| `Normal_B_cells_extended.csv` | `NBext` | Extended normal B-cell reference |
| `OncoboxCancer.csv` | `Oncobox` | Oncobox cancer reference panel |

### 3.2 Full 18-Method Table

Method names follow the `shambhala_{p_key}_{q_key}` pattern. These are the `method` slot values used in S3 keys.

| Method key | P file | Q file |
|---|---|---|
| `shambhala_P0std_Q0std` | `P0_standard.csv` | `Q0_standard.csv` |
| `shambhala_ANTE_Q0std` | `ANTE.csv` | `Q0_standard.csv` |
| `shambhala_GTExAffy_Q0std` | `GTExAffymetrix.csv` | `Q0_standard.csv` |
| `shambhala_NBBags_Q0std` | `Normal_B_BG_BAGS.csv` | `Q0_standard.csv` |
| `shambhala_NBGPL570_Q0std` | `Normal_B_GPL570.csv` | `Q0_standard.csv` |
| `shambhala_NBKass_Q0std` | `Normal_B_Kassandra.csv` | `Q0_standard.csv` |
| `shambhala_NBRNAseq_Q0std` | `Normal_B_RNASeq.csv` | `Q0_standard.csv` |
| `shambhala_NBext_Q0std` | `Normal_B_cells_extended.csv` | `Q0_standard.csv` |
| `shambhala_Oncobox_Q0std` | `OncoboxCancer.csv` | `Q0_standard.csv` |
| `shambhala_P0std_QNBKass` | `P0_standard.csv` | `Normal_B_Kassandra.csv` |
| `shambhala_ANTE_QNBKass` | `ANTE.csv` | `Normal_B_Kassandra.csv` |
| `shambhala_GTExAffy_QNBKass` | `GTExAffymetrix.csv` | `Normal_B_Kassandra.csv` |
| `shambhala_NBBags_QNBKass` | `Normal_B_BG_BAGS.csv` | `Normal_B_Kassandra.csv` |
| `shambhala_NBGPL570_QNBKass` | `Normal_B_GPL570.csv` | `Normal_B_Kassandra.csv` |
| `shambhala_NBKass_QNBKass` | `Normal_B_Kassandra.csv` | `Normal_B_Kassandra.csv` |
| `shambhala_NBRNAseq_QNBKass` | `Normal_B_RNASeq.csv` | `Normal_B_Kassandra.csv` |
| `shambhala_NBext_QNBKass` | `Normal_B_cells_extended.csv` | `Normal_B_Kassandra.csv` |
| `shambhala_Oncobox_QNBKass` | `OncoboxCancer.csv` | `Normal_B_Kassandra.csv` |

**Note on `shambhala_NBKass_QNBKass`:** P and Q are the same file (`Normal_B_Kassandra.csv`). P anchors quantile normalization shape and Q defines the target mean/std. Using the same file creates a circular normalization target but is numerically valid. Document this in the README.

---

## 4. S3 Key Conventions

### 4.1 Input: Prepared datasets

Read from the **existing** main benchmark prepared datasets — no new preparation stage is needed:

```
s3://your-bucket/FL_batch_correction/prepared/
    S0_no_removal__{imp}__exp.tsv.gz      ← full harmonizable expression matrix
    S0_no_removal__{imp}__ann.tsv.gz      ← full annotation
    {strat}__{imp}__ann.tsv.gz            ← per-strategy filtered annotation
                                            (provides sample IDs for post-hoc filtering)
```

Valid `{imp}` values: `strict`, `knn`, `softimpute`.  
Valid `{strat}` values: all 11 strategies (confirmed present in S3 for all 3 imputations).

### 4.2 Input: Calibration datasets

```
s3://your-bucket/FL_batch_correction/calibration_datasets/
    P0_standard.csv
    ANTE.csv
    GTExAffymetrix.csv
    Normal_B_BG_BAGS.csv
    Normal_B_GPL570.csv
    Normal_B_Kassandra.csv
    Normal_B_RNASeq.csv
    Normal_B_cells_extended.csv
    OncoboxCancer.csv
    Q0_standard.csv
```

All files: FL convention (rows = samples, columns = HGNC gene symbols). Load directly with `pd.read_csv(..., index_col=0)`.

### 4.3 Output: Harmonized expression

Output keys follow **exactly** the same pattern as the main framework so that the existing metrics analysis pipeline (`harmonization-metrics/`) ingests Shambhala outputs without modification:

```
s3://your-bucket/FL_batch_correction/exp/
    {strat}__{imp}__{method}__post0.tsv.gz    ← post_rm=False
    {strat}__{imp}__{method}__post1.tsv.gz    ← post_rm=True
```

Example:
```
A_confirmed_bad__knn__shambhala_P0std_Q0std__post0.tsv.gz
A_confirmed_bad__knn__shambhala_P0std_Q0std__post1.tsv.gz
```

### 4.4 Output: Metrics sidecar

Each worker writes a JSON sidecar with 22 rows (one per `strat × post_rm`):

```
s3://your-bucket/FL_batch_correction/shambhala_metrics/
    {imp}__{method}.json
```

Row keys: `strat`, `imp`, `method`, `post_rm`, `r2_batch`, `r2_diag`, `n_samples`, `n_genes`, `status`.

### 4.5 Output: Failed-jobs log (mirrored to S3)

```
s3://your-bucket/FL_batch_correction/failed_shambhala_{pod_name}.txt
```

---

## 5. New Directory Structure

```
Shambhala_containerized/
├── (existing files — unchanged)
│   ├── run_shambhala.py
│   ├── shambhala/
│   ├── octave/
│   ├── tests/
│   ├── pyproject.toml
│   └── ...
│
└── harmonization_scripts/               ← NEW
    ├── shambhala_bench_shared.py         ← shared constants, S3 utilities, strategy map
    ├── run_shambhala_job.py              ← worker: one (imp × variant) → 22 S3 outputs
    ├── run_shambhala_parallel.py         ← dispatcher: enumerate 54 jobs, launch workers
    ├── failed_jobs_shambhala.txt         ← auto-generated failure log (gitignored)
    ├── README.md                         ← cross-product usage guide
    └── k8s/
        ├── shambhala-pod.yaml            ← K8s pod manifest
        └── README.md                     ← pod launch and SSH instructions
```

---

## 6. Module Specifications

### 6.1 `shambhala_bench_shared.py`

**Purpose:** Single source of truth for all constants, S3 helpers, strategy definitions, post-removal logic, and metrics. Must NOT import `rpy2` (Shambhala is pure-Python + Octave).

**Key contents:**

```python
# ── Constants ──────────────────────────────────────────────────────────────────
S3_BUCKET       = "your-bucket"
S3_PREFIX       = "FL_batch_correction"
S3_CALIB_PREFIX = "FL_batch_correction/calibration_datasets"
BATCH_COL       = "RNA_BATCH"
BIO_COL         = "Diagnosis_cell_type_unified"

# ── Calibration dataset registry ───────────────────────────────────────────────
# Maps short key → S3 filename (relative to calibration_datasets/)
P_DATASETS: dict[str, str] = {
    "P0std":    "P0_standard.csv",
    "ANTE":     "ANTE.csv",
    "GTExAffy": "GTExAffymetrix.csv",
    "NBBags":   "Normal_B_BG_BAGS.csv",
    "NBGPL570": "Normal_B_GPL570.csv",
    "NBKass":   "Normal_B_Kassandra.csv",
    "NBRNAseq": "Normal_B_RNASeq.csv",
    "NBext":    "Normal_B_cells_extended.csv",
    "Oncobox":  "OncoboxCancer.csv",
}

Q_DATASETS: dict[str, str] = {
    "Q0std":   "Q0_standard.csv",
    "QNBKass": "Normal_B_Kassandra.csv",
}

# ── Shambhala variant registry ─────────────────────────────────────────────────
# Maps method_key → (p_short, q_short)
SHAMBHALA_VARIANTS: dict[str, tuple[str, str]] = {
    f"shambhala_{p}_{q}": (p, q)
    for q in Q_DATASETS
    for p in P_DATASETS
}  # 18 entries

# ── Strategy list (same 11 as main framework) ──────────────────────────────────
ALL_STRATEGIES: list[str] = [
    "S0_no_removal", "A_confirmed_bad", "B_extended_bad",
    "C_rnaseq_only", "D_malignant_only",
    "E1_iterative_r1", "E2_iterative_r2", "E3_iterative_r3",
    "F_microarray_only", "G_affymetrix_only", "H_affymetrix_extended",
]
ALL_IMPUTATION: list[str] = ["strict", "knn", "softimpute"]
ALL_METHODS:    list[str] = sorted(SHAMBHALA_VARIANTS.keys())

# ── S3 key functions ────────────────────────────────────────────────────────────
def s3_key_prepared_exp(strat: str, imp: str) -> str:
    return f"{S3_PREFIX}/prepared/{strat}__{imp}__exp.tsv.gz"

def s3_key_prepared_ann(strat: str, imp: str) -> str:
    return f"{S3_PREFIX}/prepared/{strat}__{imp}__ann.tsv.gz"

def s3_key_exp(strat: str, imp: str, method: str, post_rm: bool) -> str:
    suffix = "post1" if post_rm else "post0"
    return f"{S3_PREFIX}/exp/{strat}__{imp}__{method}__{suffix}.tsv.gz"

def s3_key_sidecar(imp: str, method: str) -> str:
    return f"{S3_PREFIX}/shambhala_metrics/{imp}__{method}.json"

def s3_key_calib(filename: str) -> str:
    return f"{S3_CALIB_PREFIX}/{filename}"

# ── S3 helpers ──────────────────────────────────────────────────────────────────
def s3_exists(s3_client, key: str) -> bool: ...
def download_exp_from_s3(s3_client, key: str) -> pd.DataFrame: ...
def download_ann_from_s3(s3_client, key: str) -> pd.DataFrame: ...
def upload_exp_to_s3(exp_df: pd.DataFrame, s3_client, key: str) -> None: ...

def download_calib_df(s3_client, filename: str) -> pd.DataFrame:
    # Reads FL convention CSV (rows=samples, cols=HGNC genes) from calibration_datasets/
    # Confirmed format — no transposition needed.
    ...

# ── Post-removal ────────────────────────────────────────────────────────────────
def identify_outlier_batches(
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    batch_col: str = BATCH_COL,
    n_pcs: int = 2,
    n_outliers: int = 1,
    min_batch_size: int = 20,
) -> list[str]: ...

# ── PCA R² metric ───────────────────────────────────────────────────────────────
def r2_batch(exp_df: pd.DataFrame, ann_df: pd.DataFrame,
             batch_col: str = BATCH_COL, n_pcs: int = 10) -> float: ...

# ── Memory helpers ──────────────────────────────────────────────────────────────
def free_gb() -> float: ...
```

**Implementation note:** `identify_outlier_batches`, `r2_batch`, `pca_variance_explained_by_batch`, and `free_gb` are copied verbatim from `bench_shared.py` / `run_norm_job.py`. Do **not** import from `bench_shared.py` — that module imports rpy2 at the top level, requiring a full R installation.

---

### 6.2 `run_shambhala_job.py`

**Purpose:** Worker subprocess. Given one `(imp, method)` pair, runs Shambhala on the full S0 dataset, then writes all 22 `(strat × post_rm)` outputs to S3.

**CLI signature:**
```
python run_shambhala_job.py \
    --imp         knn \
    --method      shambhala_P0std_Q0std \
    [--skip-if-exists] \
    [--out-json   /tmp/shambhala_jobs/knn__shambhala_P0std_Q0std.json] \
    [--memory-limit-gb 8.0] \
    [--n-shambhala-workers 5] \
    [--octave-bin octave] \
    [--random-seed 42] \
    [--shambhala-timeout-s 7200] \
    [--na-strategy drop]
```

**Full execution flow:**

```
Step 0 — Early-exit check (if --skip-if-exists):
    For each of 22 output keys: call s3_exists()
    If ALL 22 exist → write cached JSON rows, exit 0

Step 1 — Memory pre-flight:
    Assert free_gb() >= memory_limit_gb; exit 2 if insufficient

Step 2 — Download S0 prepared data:
    exp_full = download_exp_from_s3("S0_no_removal__{imp}__exp.tsv.gz")
    ann_full = download_ann_from_s3("S0_no_removal__{imp}__ann.tsv.gz")
    Log: shape, elapsed

Step 3 — Download calibration P and Q:
    p_key, q_key = SHAMBHALA_VARIANTS[method]
    p_df = download_calib_df(P_DATASETS[p_key])   # FL convention: rows=samples, cols=HGNC genes
    q_df = download_calib_df(Q_DATASETS[q_key])
    Log: shapes

Step 4 — Run Shambhala harmonization on full S0 expression:
    Use the shambhala Python library API directly (not subprocess):
        - Gene intersection: input ∩ P ∩ Q
        - NA handling: --na-strategy drop (default; data is already imputed — do not double-impute)
        - Gene intersection + NA strategy applied via shambhala.na_handling
        - Q statistics via shambhala.q_rescale.compute_q_statistics
        - Dispatch via shambhala.parallel.harmonize_parallel(n_workers=n_shambhala_workers)
    Log: shape, elapsed, genes_common, genes_dropped_from_{input,P,Q}
    exp_harmonized = result (FL convention: samples × genes)

    If harmonization fails: write all 22 rows with status="failed", exit 1

Step 5 — Free raw S0 expression:
    del exp_full; gc.collect()
    Log: free RAM after deletion

Step 6 — Download all 11 strategy annotations from S3:
    For each strat in ALL_STRATEGIES:
        ann_strat[strat] = download_ann_from_s3("{strat}__{imp}__ann.tsv.gz")
    Log: sample counts per strategy

Step 7 — For each (strat, post_rm) combination (22 total):
    a. Filter:
           exp_strat = exp_harmonized.loc[
               exp_harmonized.index.intersection(ann_strat[strat].index)]
           ann_s = ann_strat[strat].loc[exp_strat.index]
           Log: n_samples after filtering

    b. If post_rm=True:
           outliers = identify_outlier_batches(exp_strat, ann_s)
           if outliers:
               ann_s = ann_s[~ann_s[BATCH_COL].isin(outliers)]
               exp_strat = exp_strat.loc[ann_s.index]
           Log: outliers found (if any)

    c. Check S3 existence (if --skip-if-exists):
           key = s3_key_exp(strat, imp, method, post_rm)
           if s3_exists(key): record "cached", del exp_strat, ann_s; continue

    d. Compute metrics:
           r2_b = r2_batch(exp_strat, ann_s, BATCH_COL)
           r2_d = r2_batch(exp_strat, ann_s, BIO_COL)

    e. Upload:
           upload_exp_to_s3(exp_strat, s3, key)
           If upload fails: status="upload_failed", continue (do not abort remaining)

    f. Append row to result list:
           {"strat": strat, "imp": imp, "method": method,
            "post_rm": post_rm, "r2_batch": r2_b, "r2_diag": r2_d,
            "n_samples": len(ann_s), "n_genes": exp_strat.shape[1],
            "status": "ok"|"upload_failed"|"cached"}

    g. del exp_strat, ann_s; gc.collect()

Step 8 — Write and upload output:
    Upload JSON sidecar (22 rows) to s3_key_sidecar(imp, method)
    If --out-json: also write to local path
    Print summary: statuses, total elapsed

Step 9 — Exit 0
```

**Memory management:**
- The full harmonized expression (~200 MB–4 GB depending on imputation) stays in memory from Step 4 through all of Step 7. Step 7 only creates small filtered copies, deleted after each upload.
- `exp_full` is deleted in Step 5 before downloading strategy annotations. Never hold raw S0 expression and harmonized expression simultaneously.
- Memory budget:
  - `strict`: ~200 MB expression → 4 GB limit is safe
  - `knn`: ~2–4 GB expression → 8 GB limit required; delete `exp_full` promptly in Step 5
  - `softimpute`: similar to knn

**NA strategy for pre-imputed data:** Always pass `--na-strategy drop` (the default). The S0 knn and softimpute datasets are already imputed, so using `--na-strategy knn` would double-impute. The `drop` strategy removes any residual NaN genes before Octave and marks them NaN in the output — this is the correct behavior.

---

### 6.3 `run_shambhala_parallel.py`

**Purpose:** Dispatcher. Enumerates all 54 `(imp × method)` combinations, checks for already-finished jobs, launches `run_shambhala_job.py` subprocesses, collects results, writes failure log.

**CLI signature:**
```
python run_shambhala_parallel.py \
    [--n-workers 3] \
    [--skip-if-exists] \
    [--retry-failed] \
    [--imps strict,knn,softimpute] \
    [--methods shambhala_P0std_Q0std,shambhala_ANTE_Q0std] \
    [--n-shambhala-workers 5] \
    [--octave-bin octave] \
    [--memory-limit-gb 8.0] \
    [--timeout-s 21600] \
    [--tmp-dir /tmp/shambhala_jobs] \
    [--random-seed 42]
```

**Key design decisions:**

- **`--n-workers`**: Number of parallel `run_shambhala_job.py` subprocesses. Default: `3`. Each job runs up to 5 internal Octave workers, so 3 × 5 = 15 Octave processes maximum. Increase carefully based on pod RAM.
- **`--skip-if-exists`**: At startup, perform a single S3 `list_objects_v2` to find all existing outputs. Jobs where all 22 outputs already exist are recorded as "cached" and skipped without launching a subprocess.
- **`--timeout-s`**: Per-job timeout. Default: `21600` (6 hours). Shambhala on the full knn S0 dataset (~7,238 samples) with 5 internal workers takes approximately 2–4 hours. Set conservatively.
- **`--retry-failed`**: Re-attempts combinations listed in `failed_jobs_shambhala.txt`.
- **Memory guard**: Before launching each new subprocess, wait until `free_gb() >= memory_limit_gb`. Workers receive `--memory-limit-gb memory_limit_gb / 2`.
- **Failed-jobs log**: Local `failed_jobs_shambhala.txt` + S3 mirror at `FL_batch_correction/failed_shambhala_{POD_NAME}.txt`. Format: one `{imp}__{method}` key per line.
- **Metrics upload**: At completion, reads all per-job JSON sidecars and uploads consolidated `shambhala_metrics.csv` to `s3://your-bucket/FL_batch_correction/shambhala_metrics.csv`.

**Internal flow (mirroring `run_norm_parallel.py`):**

```
1. Sync failed-jobs log from S3 (merge with local)
2. List all S3 outputs with list_objects_v2 to find what already exists
3. Build job list: all (imp, method) pairs − already_complete − skip_keys
4. Print: N jobs, M workers, timeout, memory guard settings
5. ThreadPoolExecutor(max_workers=N):
       for each (imp, method):
           wait_for_memory(memory_limit_gb)
           launch run_shambhala_job.py subprocess
           stream stdout with prefix [imp × method]
           wait up to timeout_s; kill on timeout and record as failed
           read JSON sidecar → collect rows
6. Update and save failed-jobs log
7. Sync failed-jobs log to S3
8. Read all per-job JSON sidecars → build DataFrame → upload shambhala_metrics.csv
```

**Exit codes from `run_shambhala_job.py`:**

| Code | Meaning |
|---|---|
| 0 | Success (or all outputs already cached) |
| 1 | Shambhala harmonization failed |
| 2 | Insufficient RAM (pre-flight check failed) |
| 3 | Prepared dataset not found in S3 |

---

## 7. Data Flow Diagram

```
S3: FL_batch_correction/prepared/
    S0_no_removal__{imp}__exp.tsv.gz        ──┐
    S0_no_removal__{imp}__ann.tsv.gz        ──┤
    {strat}__{imp}__ann.tsv.gz (×11)        ──┤
                                              │
S3: FL_batch_correction/calibration_datasets/ │
    {P_FILE}.csv                            ──┤
    {Q_FILE}.csv                            ──┘
                                              │
                         run_shambhala_job.py │
                         ┌────────────────────┴──────────────────────┐
                         │                                           │
                         │  1. Download S0 expression               │
                         │  2. Download P, Q                        │
                         │  3. Shambhala harmonize (ALL samples)    │
                         │  4. del exp_full → free RAM              │
                         │  5. Download 11 strategy annotations     │
                         │  6. For each (strat × post_rm):          │
                         │      filter → [outlier removal] → upload │
                         │  7. Write metrics sidecar                │
                         └──────────────────────────────────────────┘
                                              │
S3: FL_batch_correction/exp/
    {strat}__{imp}__shambhala_{p}_{q}__post{0|1}.tsv.gz   (22 files per job)

S3: FL_batch_correction/shambhala_metrics/
    {imp}__{method}.json                    (per-job sidecar)
    shambhala_metrics.csv                   (consolidated — written by dispatcher)
```

---

## 8. K8s Pod Specification (`k8s/shambhala-pod.yaml`)

### 8.1 Key differences from `harmonization-scripts/k8s/pod-ssh.yaml`

| Aspect | Main benchmark pod | Shambhala pod |
|---|---|---|
| Base image | `ubuntu:24.04` | `ubuntu:24.04` (same) |
| R requirement | R 4.5.x + Bioconductor (~30 min) | **None** |
| Python version | 3.11 (for reComBat compat) | **3.11** (same, for consistency) |
| Key Python packages | rpy2, harmonypy, scanorama, inmoose, reComBat… | boto3, pandas, numpy, scikit-learn, psutil (only) |
| Octave | Required (for Shambhala) | **Required** (same Octave + matlab wrapper) |
| Startup time | ~30–40 min (R packages) | **~5–10 min** (no R) |
| Pod name | `shambhala-pvc` | `shambhala-bench` |
| CPU/RAM | 16 vCPU / 128 GiB | **16 vCPU / 128 GiB** |
| Working directory | `/app/harmonization-scripts` | `/app/Shambhala_containerized` |
| PVC | `shambhala-pvc` | `shambhala-pvc` (same — shared `/workspace`) |

### 8.2 ConfigMap: `shambhala-deps` (requirements.txt)

```yaml
data:
  requirements.txt: |
    boto3>=1.28.17
    botocore>=1.31.17
    pandas>=1.5.3,<2.0
    numpy>=1.24,<2.0
    scikit-learn>=1.3
    psutil>=5.9
    awscli>=1.29
```

**No R, rpy2, harmonypy, scanorama, inmoose, or reComBat.** The Shambhala pipeline is pure Python + Octave.

### 8.3 Startup sequence (pod args)

```bash
# -- 1. SSH (done first so you can connect and monitor) --
apt-get update -qq
apt-get install -y -qq openssh-server rsync curl wget gnupg ca-certificates tmux htop
# [standard SSH setup identical to main pod]
/usr/sbin/sshd -D &
echo "[startup] SSH ready"

# -- 2. System libs for Python and Octave --
apt-get install -y -qq build-essential python3 python3-pip python3-dev python-is-python3 \
    octave octave-statistics liboctave-dev
# matlab shim (Shambhala2 calls system("matlab ...") — strip MATLAB-only flags)
printf '#!/bin/bash\nargs=()\nfor a in "$@"; do\n  case "$a" in\n    -nodesktop|-nosplash|-nodisplay) ;;\n    *) args+=("$a") ;;\n  esac\ndone\nexec octave --no-gui "${args[@]}"\n' \
    > /usr/local/bin/matlab && chmod +x /usr/local/bin/matlab
echo "[startup] System libs + Octave installed"

# -- 3. Python packages (fast — no R) --
pip3 install --break-system-packages --no-cache-dir -r /etc/shambhala-deps/requirements.txt
echo "[startup] Python packages installed"

# -- 4. Working directory --
mkdir -p /app/Shambhala_containerized
echo 'cd /app/Shambhala_containerized' >> /root/.bashrc
echo "[startup] Environment ready. Sync code with rsync and run."
sleep infinity
```

**No R installation step.** Startup completes in ~5–10 minutes (vs 30–40 for the main pod).

### 8.4 Resource spec

```yaml
resources:
  requests:
    cpu: "16"
    memory: "128Gi"
  limits:
    cpu: "16"
    memory: "128Gi"
```

**Justification:** With `--n-workers 3` and `--n-shambhala-workers 5`, peak RAM = 3 jobs × ~4 GB (knn expression) + 15 Octave processes × ~1 GB = ~27 GB. 128 GiB provides comfortable headroom for dispatcher overhead and simultaneous S3 uploads.

---

## 9. `README.md` for `harmonization_scripts/`

The README must cover:

1. **What this does** — cross-product description, 54 runs → 1,188 outputs.
2. **Prerequisites** — pod running, code synced via rsync, S3 credentials mounted.
3. **Quick start** — copy-paste commands for common runs.
4. **Full command reference** — all CLI flags for `run_shambhala_parallel.py` and `run_shambhala_job.py`.
5. **Method naming reference** — the 18 variant table (P key × Q key → method name).
6. **Resuming a partial run** — `--skip-if-exists`, `--retry-failed`, failure log location and format.
7. **Monitoring** — `tail -f /workspace/shambhala_run.log`, estimating progress from S3 object count.
8. **Timing estimates** — per-variant runtime table:
   - strict: ~30–60 min (small gene set, few NAs)
   - knn: ~2–4 h (large gene set, all samples)
   - softimpute: ~2–4 h (similar to knn)
9. **Output format** — same `.tsv.gz` format as main benchmark; compatible with `load_cross_product_results.py`.
10. **NA strategy note** — always use `--na-strategy drop` (default) since S0 data is pre-imputed.

Key command examples:
```bash
# Full run (all 54 jobs)
nohup python harmonization_scripts/run_shambhala_parallel.py \
    --n-workers 3 \
    --n-shambhala-workers 5 \
    --octave-bin /usr/bin/octave \
    --skip-if-exists \
    --memory-limit-gb 10.0 \
    --timeout-s 21600 \
    --random-seed 42 \
    > /workspace/shambhala_run.log 2>&1 &

# Smoke test: 1 variant, strict imputation only
python harmonization_scripts/run_shambhala_job.py \
    --imp strict \
    --method shambhala_P0std_Q0std \
    --n-shambhala-workers 2 \
    --octave-bin /usr/bin/octave \
    --random-seed 42 \
    --out-json /tmp/test_result.json

# Targeted subset (e.g., 2 variants, strict only)
python harmonization_scripts/run_shambhala_parallel.py \
    --imps strict \
    --methods shambhala_P0std_Q0std,shambhala_NBBags_Q0std \
    --n-workers 2 \
    --n-shambhala-workers 3 \
    --skip-if-exists \
    --octave-bin /usr/bin/octave
```

---

## 10. `k8s/README.md` for Pod Operations

The K8s README must cover step by step:

1. **Prerequisites** — `kubectl` configured for the cluster, `~/.aws/credentials` present, PVC `shambhala-pvc` exists.

2. **One-time setup** (if not already done):
   ```bash
   kubectl create secret generic aws-credentials \
       --from-file=credentials=$HOME/.aws/credentials \
       --from-file=config=$HOME/.aws/config \
       --namespace <your-namespace>
   ```

3. **Deploy the pod:**
   ```bash
   kubectl apply -f harmonization_scripts/k8s/shambhala-pod.yaml -n <your-namespace>
   ```

4. **Verify pod creation and monitor startup:**
   ```bash
   kubectl get pod shambhala-bench -n <your-namespace>   # check Pending → Running
   kubectl logs -f shambhala-bench -n <your-namespace>   # wait for "Environment ready."
   ```
   Startup takes ~5–10 min. SSH becomes available within ~1 min (before Python/Octave finish installing).

5. **Port-forward and SSH access:**
   ```bash
   kubectl port-forward pod/shambhala-bench 2223:22 -n <your-namespace> &
   ssh -p 2223 -i ~/.ssh/id_ed25519 -o StrictHostKeyChecking=no root@localhost
   ```
   Or add to `~/.ssh/config`:
   ```
   Host shambhala-pod
       HostName localhost
       Port 2223
       User root
       IdentityFile ~/.ssh/id_ed25519
       StrictHostKeyChecking no
   ```
   Then: `ssh shambhala-pod`

6. **Sync code from local Mac to pod:**
   ```bash
   rsync -avz --progress \
       -e "ssh -p 2223 -i ~/.ssh/id_ed25519 -o StrictHostKeyChecking=no" \
       /Users/user890/Desktop/fl_subset/shambhala_adoption/Shambhala_containerized/ \
       root@localhost:/app/Shambhala_containerized/
   ```

7. **Install Shambhala package inside pod** (after sync):
   ```bash
   cd /app/Shambhala_containerized
   pip3 install --break-system-packages -e .
   ```

8. **Run the benchmark** — see `harmonization_scripts/README.md` for full command reference.

9. **Check progress:**
   ```bash
   tail -f /workspace/shambhala_run.log     # live log
   aws s3 ls s3://your-bucket/FL_batch_correction/exp/ \
       | grep shambhala | wc -l             # S3 output count (target: 1188)
   htop                                     # CPU / RAM usage
   ```

10. **Teardown:**
    ```bash
    kubectl delete pod shambhala-bench -n <your-namespace>
    ```
    The PVC `/workspace` and S3 outputs persist after pod deletion.

---

## 11. Integration with Main Benchmark Framework

Because Shambhala outputs use **identical S3 key patterns** as the main benchmark, they are automatically compatible with:

- **`load_cross_product_results.py`** — `load_metrics()`, `load_exp()`, `load_top_k()` work without modification.
- **`harmonization-metrics/`** — the metrics pipeline consumes any key matching `exp/{strat}__{imp}__{method}__post{0|1}.tsv.gz`. Shambhala outputs will be picked up in the next metrics run without any changes.
- **`harmonization_benchmark.ipynb`** — PCA R² heatmaps and ranking plots will include all 18 Shambhala variants automatically.

**The only optional change in the main framework:** add Shambhala method keys to `ALL_METHODS` in `run_norm_parallel.py` if you ever want the main dispatcher to check/skip them. Not required — Shambhala is run by its own dispatcher.

---

## 12. Known Limitations and Edge Cases

| Issue | Description | Handling |
|---|---|---|
| Gene coverage mismatch | P or Q may cover a different gene set than the S0 expression; intersection is taken at Step 4 | Logged: genes_common, genes_dropped_from_{input,P,Q} |
| `shambhala_NBKass_QNBKass` | P == Q (same file); circular normalization target | Numerically valid; documented in README |
| knn/softimpute imputation size | S0 knn/softimpute expression may be 2–4 GB uncompressed | `del exp_full; gc.collect()` in Step 5; `--memory-limit-gb 10` recommended |
| Octave binary path | `octave` may not be on PATH inside pod | Pass `--octave-bin /usr/bin/octave` explicitly in all pod commands |
| Shambhala NA strategy vs pre-imputed input | knn/softimpute data is already imputed — do not double-impute | Default `--na-strategy drop` is correct; documented in README |
| Parallel job RAM competition | `--n-workers 3` × knn jobs = ~12 GB peak simultaneously | Dispatcher memory guard enforces minimum free RAM before each launch |
| Octave k-means non-reproducibility | k-means initialization uses random state | Always pass `--random-seed 42` for reproducible benchmark |
| Existing Shambhala bugs (unfixed) | 3 bugs in `run_shambhala.py` / `octave_bridge.py` affect 8/30 tests | **Must be fixed before implementation** — see §13 Phase 0 |

---

## 13. Implementation Phases Overview

See §14 for the detailed task-level TODO list.

| Phase | What | Output |
|---|---|---|
| 0 | Fix 3 existing bugs in `run_shambhala.py` + `octave_bridge.py` | All 30 tests pass |
| 1 | `shambhala_bench_shared.py` | Core library, unit-testable without S3 |
| 2 | `run_shambhala_job.py` | Worker, testable on local data |
| 3 | `run_shambhala_parallel.py` | Dispatcher, mirrors `run_norm_parallel.py` |
| 4 | `k8s/shambhala-pod.yaml` + `k8s/README.md` | Pod manifest and ops guide |
| 5 | `harmonization_scripts/README.md` | User-facing usage guide |
| 6 | End-to-end validation on pod | 1,188 S3 outputs confirmed |

---

## 14. Detailed TODO List

### Phase 0 — Fix Existing Shambhala Bugs (prerequisite for everything)

These three bugs are documented in `CLAUDE.md` as unfixed. **8/30 tests currently fail.** All three must be fixed before cross-product scripts can rely on the library.

- [x] **0.1** Fix Bug 1 in `shambhala/octave_bridge.py:203`:
  - Change `f"source('Shambhala2_piped.m');"` → `f"source('{octave_scripts_dir}/Shambhala2_piped.m');"`
  - This resolves the "no such file" Octave error affecting all 8 failing tests.
  - **Status:** Already fixed in a prior session.

- [x] **0.2** Fix Bug 2 in `run_shambhala.py` Step 5:
  - `apply_na_strategy` must receive FL convention (samples×genes), not genes×samples.
  - Change Step 5 to: `clean_df, na_mask = na_handling.apply_na_strategy(input_T.T, ...)` then `clean_T = clean_df.T`.
  - Align P: `p_T = p_T.loc[clean_T.index]`.
  - **Status:** Already fixed in a prior session.

- [x] **0.3** Fix Bug 3 in `run_shambhala.py` Step 6:
  - `compute_q_statistics` must receive FL convention (samples×genes).
  - Change Step 6 to: `q_df_filtered = q_T.loc[clean_T.index].T` then `rm, rs = q_rescale.compute_q_statistics(q_df_filtered, ...)`.
  - **Status:** Already fixed in a prior session.

- [x] **0.4** Fix same Bug 2 + Bug 3 pattern in `tests/test_integration.py::_run_pipeline`.
  - **Status:** Already fixed in a prior session.

- [x] **0.5** Fix Bug 2 in `run_shambhala.py` Step 8:
  - `restore_na_genes` must receive FL convention.
  - Change Step 8 to: `result_df = na_handling.restore_na_genes(harmonized_T.T, na_mask)`.
  - **Status:** Already fixed in a prior session.

- [x] **0.6** Run full test suite and confirm all 30 tests pass:
  ```bash
  cd /home/jovyan/Projects/Follicular_lymphoma_disser/shambhala_adoption/Shambhala_containerized
  pytest tests/ -v
  ```
  Expected: 30 passed, 0 failed (some may still skip if Octave is unavailable in the test environment).
  **Result: 30 passed in 659s — all tests pass.**

- [x] **0.7** Run the CLI smoke test end-to-end to confirm the fixed pipeline works:
  ```bash
  python run_shambhala.py \
      --input tests/fixtures/input_10samples.csv \
      --P tests/fixtures/P0_small.csv \
      --Q tests/fixtures/Q0_small.csv \
      --output /tmp/harmonized_test.csv \
      --octave-bin /opt/conda/bin/octave \
      --n-workers 1 --random-seed 42
  ```
  **Status:** Covered by test_integration.py::test_single_worker_matches_reference — passes.

---

### Phase 1 — Core Library `shambhala_bench_shared.py`

- [x] **1.1** Create `harmonization_scripts/` directory.

- [x] **1.2** Create `harmonization_scripts/shambhala_bench_shared.py` with all constants:
  - `S3_BUCKET`, `S3_PREFIX`, `S3_CALIB_PREFIX`, `BATCH_COL`, `BIO_COL`
  - `P_DATASETS` dict (9 entries)
  - `Q_DATASETS` dict (2 entries)
  - `SHAMBHALA_VARIANTS` dict comprehension (18 entries)
  - `ALL_STRATEGIES` list (11 entries — identical to main framework)
  - `ALL_IMPUTATION` list (`["strict", "knn", "softimpute"]`)
  - `ALL_METHODS` list (sorted SHAMBHALA_VARIANTS keys)

- [x] **1.3** Add S3 key functions:
  - `s3_key_prepared_exp(strat, imp) -> str`
  - `s3_key_prepared_ann(strat, imp) -> str`
  - `s3_key_exp(strat, imp, method, post_rm) -> str`
  - `s3_key_sidecar(imp, method) -> str`
  - `s3_key_calib(filename) -> str`

- [x] **1.4** Add S3 I/O helpers (copy + adapt from `bench_shared.py`):
  - `s3_exists(s3_client, key) -> bool`
  - `download_exp_from_s3(s3_client, key) -> pd.DataFrame` (handles `.tsv.gz`)
  - `download_ann_from_s3(s3_client, key) -> pd.DataFrame` (handles `.tsv.gz`)
  - `upload_exp_to_s3(exp_df, s3_client, key) -> None` (gzip TSV)
  - `download_calib_df(s3_client, filename) -> pd.DataFrame` (CSV, FL convention, no transpose)

- [x] **1.5** Add metric and post-removal functions (copy verbatim from `bench_shared.py` — no rpy2):
  - `pca_variance_explained_by_batch(exp_df, batch_series, n_components=10) -> float`
  - `r2_batch(exp_df, ann_df, batch_col=BATCH_COL, n_pcs=10) -> float`
  - `identify_outlier_batches(exp_df, ann_df, ...) -> list[str]`

- [x] **1.6** Add `free_gb() -> float` (copy from `run_norm_job.py`).

- [x] **1.7** Verify: `python -c "from shambhala_bench_shared import SHAMBHALA_VARIANTS, ALL_STRATEGIES; assert len(SHAMBHALA_VARIANTS) == 18; assert len(ALL_STRATEGIES) == 11; print('OK')"` → `OK`.

- [x] **1.8** Verify S3 key format: confirm one generated key matches the expected pattern:
  ```python
  from shambhala_bench_shared import s3_key_exp
  assert s3_key_exp("A_confirmed_bad", "knn", "shambhala_P0std_Q0std", False) == \
      "FL_batch_correction/exp/A_confirmed_bad__knn__shambhala_P0std_Q0std__post0.tsv.gz"
  ```

---

### Phase 2 — Worker `run_shambhala_job.py`

- [x] **2.1** Write `harmonization_scripts/run_shambhala_job.py` argument parser:
  - `--imp` (required)
  - `--method` (required)
  - `--skip-if-exists` (flag)
  - `--out-json` (optional path)
  - `--memory-limit-gb` (float, default 8.0)
  - `--n-shambhala-workers` (int, default 5)
  - `--octave-bin` (str, default "octave")
  - `--random-seed` (int, optional)
  - `--shambhala-timeout-s` (int, default 7200)
  - `--na-strategy` (choices: drop/knn, default "drop")

- [x] **2.2** Implement Step 0: early-exit check.
  - Build list of all 22 S3 keys for this `(imp, method)`.
  - If `--skip-if-exists` and all 22 exist: write cached JSON rows, exit 0.

- [x] **2.3** Implement Step 1: memory pre-flight (`free_gb() >= memory_limit_gb`; exit 2 on failure).

- [x] **2.4** Implement Step 2: download S0 expression and annotation.
  - Use `s3_key_prepared_exp("S0_no_removal", imp)` and `s3_key_prepared_ann("S0_no_removal", imp)`.
  - Log shape and download time.

- [x] **2.5** Implement Step 3: download calibration P and Q.
  - Look up `p_key, q_key = SHAMBHALA_VARIANTS[method]`.
  - Download P and Q using `download_calib_df`.
  - Log shapes.

- [x] **2.6** Implement Step 4: run Shambhala harmonization using the Python library API.
  - Import `shambhala.na_handling`, `shambhala.q_rescale`, `shambhala.parallel`.
  - Follow the gene intersection + NA strategy + Q statistics + `harmonize_parallel` pattern from `run_shambhala.py` main().
  - Wrap in try/except: on failure write 22 `status="failed"` rows, exit 1.
  - Log: genes common, genes dropped from each source, elapsed.

- [x] **2.7** Implement Step 5: `del exp_full; gc.collect()`. Log free RAM after deletion.

- [x] **2.8** Implement Step 6: download all 11 strategy annotations.
  - Loop `ALL_STRATEGIES`, call `download_ann_from_s3(s3_key_prepared_ann(strat, imp))`.
  - Store in `ann_strat: dict[str, pd.DataFrame]`.
  - Log sample counts per strategy.

- [x] **2.9** Implement Step 7 loop over 22 `(strat, post_rm)` combinations:
  - Filter `exp_harmonized` to strategy sample IDs.
  - If `post_rm=True`: call `identify_outlier_batches`; filter.
  - Skip if `--skip-if-exists` and key already exists.
  - Compute `r2_batch` and `r2_diag`.
  - Upload; handle `upload_failed` without aborting.
  - Delete filtered copy after each upload.

- [x] **2.10** Implement Step 8: write JSON sidecar to local path + upload to S3.

- [ ] **2.11** Local smoke test (no pod required — use test fixtures):
  ```bash
  python harmonization_scripts/run_shambhala_job.py \
      --imp strict \
      --method shambhala_P0std_Q0std \
      --n-shambhala-workers 1 \
      --octave-bin /opt/conda/bin/octave \
      --random-seed 42 \
      --out-json /tmp/test_job.json
  ```
  Confirm: JSON sidecar has 22 rows; S3 keys are created.
  **Requires live S3 access — run on pod (Phase 6).**

- [ ] **2.12** Verify exit codes: simulate S3 missing (exit 3) by running with a non-existent imp name.
  **Requires live S3 access — run on pod (Phase 6).**

---

### Phase 3 — Dispatcher `run_shambhala_parallel.py`

- [x] **3.1** Write `harmonization_scripts/run_shambhala_parallel.py` argument parser (all flags from §6.3).

- [x] **3.2** Implement `_load_failed_log()` and `_save_failed_log()` using local `failed_jobs_shambhala.txt`.

- [x] **3.3** Implement `_sync_failed_log_from_s3()` and `_sync_failed_log_to_s3()` (mirrors `run_norm_parallel.py`).

- [x] **3.4** Implement S3 existence scan at startup: `list_objects_v2` on `exp/` prefix, build set of existing output keys.

- [x] **3.5** Implement job enumeration and filtering:
  - All `(imp, method)` pairs in `ALL_IMPUTATION × ALL_METHODS`.
  - Remove pairs where all 22 outputs already exist (if `--skip-if-exists`).
  - Remove pairs in `skip_keys` (failed log).

- [x] **3.6** Implement `_wait_for_memory(min_free_gb)` and `_free_gb()`.

- [x] **3.7** Implement `ThreadPoolExecutor` dispatch:
  - Launch `run_shambhala_job.py` subprocess for each job.
  - Forward stdout with `[imp × method]` prefix.
  - Handle timeout (kill + record failed).
  - Read JSON sidecar and collect rows on success.

- [x] **3.8** Implement final metrics consolidation:
  - Gather all JSON sidecars into a single DataFrame.
  - Upload `shambhala_metrics.csv` to S3.

- [x] **3.9** Wire SIGTERM handler to sync failed-jobs log to S3.

- [ ] **3.10** Test with 1 worker, 1 imp, 1 method, `--skip-if-exists`:
  ```bash
  python harmonization_scripts/run_shambhala_parallel.py \
      --n-workers 1 \
      --imps strict \
      --methods shambhala_P0std_Q0std \
      --skip-if-exists \
      --octave-bin /opt/conda/bin/octave
  ```
  **Requires live S3 access — run on pod (Phase 6).**

---

### Phase 4 — K8s Pod Manifest and README

- [x] **4.1** Create `harmonization_scripts/k8s/` directory.

- [x] **4.2** Write `harmonization_scripts/k8s/shambhala-pod.yaml`:
  - ConfigMap `shambhala-ssh-pubkey` (same SSH key as main pod).
  - ConfigMap `shambhala-deps` (requirements.txt — Python only, no R).
  - Pod `shambhala-bench` with:
    - `ubuntu:24.04` image.
    - 4-step startup: SSH → Python+Octave → pip install → working dir.
    - `matlab` shim script (same as main pod).
    - 16 vCPU / 128 GiB resources.
    - Volumes: PVC `shambhala-pvc`, aws-credentials, dshm, ssh-pubkey, shambhala-deps.
    - Tolerations and affinity for `node-group=<your-namespace>` (same as main pod).
    - `POD_NAME` env var via fieldRef.

- [x] **4.3** Write `harmonization_scripts/k8s/README.md` (full step-by-step guide as specified in §10):
  - Prerequisites.
  - One-time secret creation.
  - `kubectl apply` + `kubectl get pod` + `kubectl logs -f`.
  - Port-forward + SSH config entry.
  - rsync sync command.
  - `pip install -e .` inside pod.
  - Run benchmark pointer.
  - Progress monitoring commands.
  - Teardown.

- [ ] **4.4** Deploy pod to cluster and verify startup completes in <15 min:
  ```bash
  kubectl apply -f harmonization_scripts/k8s/shambhala-pod.yaml -n <your-namespace>
  kubectl get pod shambhala-bench -n <your-namespace>
  kubectl logs -f shambhala-bench -n <your-namespace>
  ```
  Wait for `[startup] Environment ready.`

- [ ] **4.5** SSH in and verify Octave is functional:
  ```bash
  ssh shambhala-pod
  octave --version          # should print Octave version
  python3 --version         # should print Python 3.x
  which matlab              # should print /usr/local/bin/matlab
  ```

---

### Phase 5 — `harmonization_scripts/README.md`

- [x] **5.1** Write `harmonization_scripts/README.md` covering all 10 sections listed in §9.

- [x] **5.2** Include the full 18-method variant table.

- [x] **5.3** Include timing estimates table (strict ~30–60 min, knn/softimpute ~2–4 h per variant).

- [x] **5.4** Include explicit note: always pass `--octave-bin /usr/bin/octave` inside the pod.

- [x] **5.5** Include note about NA strategy: use `--na-strategy drop` (default) for pre-imputed S0 data.

---

### Phase 6 — End-to-End Validation on Pod

- [ ] **6.1** Sync code to pod:
  ```bash
  rsync -avz --progress \
      -e "ssh -p 2223 -i ~/.ssh/id_ed25519 -o StrictHostKeyChecking=no" \
      /Users/user890/Desktop/fl_subset/shambhala_adoption/Shambhala_containerized/ \
      root@localhost:/app/Shambhala_containerized/
  ```

- [ ] **6.2** Install the package inside the pod:
  ```bash
  cd /app/Shambhala_containerized && pip3 install --break-system-packages -e .
  ```

- [ ] **6.3** Run smoke test (1 variant, strict imputation):
  ```bash
  python harmonization_scripts/run_shambhala_job.py \
      --imp strict \
      --method shambhala_P0std_Q0std \
      --n-shambhala-workers 2 \
      --octave-bin /usr/bin/octave \
      --random-seed 42 \
      --out-json /workspace/smoke_test.json
  ```

- [ ] **6.4** Verify 22 S3 keys created:
  ```bash
  aws s3 ls s3://your-bucket/FL_batch_correction/exp/ \
      | grep "shambhala_P0std_Q0std" | wc -l   # expect 22
  ```

- [ ] **6.5** Verify JSON sidecar has 22 rows and metrics are plausible (r2_batch < 0.95 for most strategies).

- [ ] **6.6** Load one output and confirm shape matches S0 strict input:
  ```python
  import pandas as pd
  df = pd.read_csv(
      "s3://your-bucket/FL_batch_correction/exp/"
      "S0_no_removal__strict__shambhala_P0std_Q0std__post0.tsv.gz",
      sep="\t", index_col=0, compression="gzip"
  )
  print(df.shape)  # should match S0 strict prepared dataset shape
  ```

- [ ] **6.7** Run full benchmark (all 54 jobs):
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

- [ ] **6.8** Monitor progress until completion:
  ```bash
  tail -f /workspace/shambhala_run.log
  aws s3 ls s3://your-bucket/FL_batch_correction/exp/ \
      | grep shambhala | wc -l   # target: 1188
  ```

- [ ] **6.9** Confirm `shambhala_metrics.csv` is uploaded to S3 and contains 1,188 rows (or 54 × 22 minus any cached/failed).

- [x] **6.10** Update `Shambhala_containerized/CLAUDE.md` "Known Bugs" section to reflect that the 3 bugs are fixed and all tests pass.
  **Status:** Done — section renamed "Bug History (all fixed as of 2026-05-15)" in prior session.

---

## 15. File Size Estimates

| Output file | Approx. size (compressed) |
|---|---|
| `S0_no_removal__strict__shambhala_P0std_Q0std__post0.tsv.gz` | ~180 MB |
| `S0_no_removal__knn__shambhala_P0std_Q0std__post0.tsv.gz` | ~500 MB |
| `A_confirmed_bad__knn__shambhala_P0std_Q0std__post0.tsv.gz` | ~450 MB |
| All 396 outputs (strict only: 18 variants × 11 × 2) | ~7 GB |
| All 396 outputs (knn only) | ~24 GB |
| All 1,188 outputs (all 3 imputations) | ~35–40 GB |

All final outputs go directly to S3; only temporary per-job JSON sidecars land on `/workspace`.
