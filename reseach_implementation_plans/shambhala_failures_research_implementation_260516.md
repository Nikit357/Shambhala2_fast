# Shambhala Pod Failures — Research & Implementation Plan

**Date:** 2026-05-16  
**Log analyzed:** `logs/shambhala_pod_logs_260516.txt`  
**Status:** Root causes identified; fixes specified below.

---

## 1. What Happened — Overview

The benchmark was launched at 05:25:40 with 54 jobs (3 imputations × 18 P/Q variants), 3 parallel workers, 5 Shambhala workers per job, 21600 s timeout.

The log has 197 lines covering the first ~5 minutes. This is not a crash of the dispatcher — the run is still ongoing. Only the first batch of 3 jobs is visible:

| Job | Outcome | Duration |
|---|---|---|
| `strict×shambhala_GTExAffy_Q0std` | FAILED (exit=1) | 25 s |
| `strict×shambhala_ANTE_Q0std` | Running | still in progress |
| `strict×shambhala_ANTE_QNBKass` | Running | still in progress |
| `strict×shambhala_GTExAffy_QNBKass` | FAILED (exit=1) | 24 s (next slot) |
| `strict×shambhala_NBBags_Q0std` | Starting | log truncated |

Both GTExAffy variants failed within 25 seconds with an Octave crash. The ANTE variants appear to have continued (they had non-zero gene intersections). The log was captured while the run was still in progress.

---

## 2. Issue 1 — GTExAffy P Dataset Has Wrong Data Orientation (CRITICAL)

### Symptom

```
Gene intersection: 0 common (dropped 3447 from input, 651 from P, 11887 from Q)
NA strategy 'drop': 0 genes dropped.
Q statistics computed for 0 genes.
```

Both `GTExAffy_Q0std` and `GTExAffy_QNBKass` hit this. The P shape log line shows:

```
P: (20254, 651)   ← as read by download_calib_df
```

### Root Cause

`download_calib_df` reads the calibration CSV with FL convention (rows = samples, columns = genes):

```python
return pd.read_csv(buf, index_col=0)  # expects: rows=samples, cols=HGNC genes
```

The `GTExAffymetrix.csv` file is stored in **genes × samples** orientation (the original Shambhala2 format, not FL convention):

- 20254 rows → GPL570 Affymetrix gene/probe IDs — but `read_csv` treats them as **sample IDs**
- 651 columns → actual GTEx sample identifiers (e.g., `GTEX-1A2B3-0006-SM-...`) — but `read_csv` treats them as **gene names**

Since the 651 "gene names" are GTEx sample identifiers, their intersection with 3447 HGNC gene symbols from S0 is exactly **zero**.

This is a data preparation error: the file was uploaded to S3 in genes × samples format without transposing to FL convention. The shape `(20254, 651)` is consistent with ~20,000 Affymetrix GPL570 genes × ~651 GTEx microarray samples.

### Impact

- All 6 GTExAffy jobs (3 imputations × 2 Q variants) will fail. This accounts for 132 out of 1,188 planned S3 outputs being permanently missed.
- The failure mode is an Octave crash (Issue 2 below), not a graceful error.

### Fix

Re-prepare `GTExAffymetrix.csv` with FL convention (samples × HGNC gene symbols) and re-upload to S3 at `FL_batch_correction/calibration_datasets/GTExAffymetrix.csv`.

Steps:
1. Download the raw file: `aws s3 cp s3://your-bucket/FL_batch_correction/calibration_datasets/GTExAffymetrix.csv ./GTExAffymetrix_raw.csv`
2. Inspect: `pd.read_csv('GTExAffymetrix_raw.csv', index_col=0).shape` — confirm it is `(20254, 651)`.
3. Transpose: `df.T` — gives (651, 20254) = 651 samples × 20254 genes.
4. Check that column names are now HGNC symbols (not GTEx sample IDs). If they are probe IDs like `201820_at`, apply the same GPL570 → HGNC mapping used in `GPL570_mapping.py` in the parent project.
5. Align to the common gene set and re-upload.

**Verification test** (run locally before re-upload): `download_calib_df(s3, 'GTExAffymetrix.csv').columns[:5]` should return recognizable HGNC gene symbols (e.g., `BRCA1`, `TP53`).

---

## 3. Issue 2 — No Guard Against Empty Gene Intersection → Octave Crash (CRITICAL)

### Symptom

The Octave crash that follows the 0-gene intersection:

```
error: r_sort(5): out of bound 0 (dimensions are 0x1)
error: called from
    kmeans at line 11 column 14
    CuBlock at line 50 column 15
    Shambhala2_piped at line 25 column 11
```

### Root Cause Trace

When `common_genes` is empty (`len == 0`), the pipeline does not abort. It continues:

```python
# run_shambhala_job.py — no check here
input_T = input_T.loc[common_genes]   # → shape (0, 7174) — zero genes, 7174 samples
p_T = p_T.loc[clean_T.index]          # → shape (0, 651) — zero genes
```

A 0-gene matrix is then serialized as a TSV with 0 data rows and piped to Octave. Inside Octave:

```matlab
% Shambhala2_piped.m line 25
dataN = CuBlock(real(EXP), [], k);    % EXP is (0 × 1+651) → 0 rows

% CuBlock.m line 50
indProbes = kmeans(data, k, 'maxiter', 1000);  % data is (0, 652)

% kmeans.m line 6, 10, 11
[n_probes, n_samples] = size(data);   % n_probes = 0
[~, r_sort] = sort(rand(0, 1));       % r_sort is 0×1 (empty)
rand_idx = r_sort(1:k);               % r_sort(1:5) → OUT OF BOUNDS on 0-element vector
```

Error message: `r_sort(5): out of bound 0 (dimensions are 0x1)` — exactly matches the trace.

### Fix A — Guard in `run_shambhala_job.py`

Add after the gene intersection block (line ~228):

```python
if len(common_genes) == 0:
    raise ValueError(
        f"Gene intersection is empty for method={method!r}. "
        f"P dataset '{P_DATASETS[p_key]}' gene names do not overlap with input S0 genes. "
        "Check that the calibration CSV uses HGNC gene symbols (FL convention: rows=samples, cols=genes)."
    )
```

This converts a cryptic Octave crash into a clear Python error with diagnostic context.

### Fix B — Guard in `run_shambhala.py` (CLI)

Same guard after step 4 (gene intersection), before step 5:

```python
if len(common_genes) == 0:
    logger.error(
        "Gene intersection is empty. P and Q datasets must share gene symbols with the input. "
        "Check that all inputs use the same gene identifier format (HGNC symbols)."
    )
    sys.exit(1)
```

### Fix C — Guard in `kmeans.m` (defensive)

Add a size check at the top of `kmeans.m` so an Octave-level crash becomes an informative Octave error:

```matlab
function indProbes = kmeans(data, k, varargin)
    [n_probes, n_samples] = size(data);
    if n_probes == 0
        error('kmeans: input data has 0 rows (0 genes). Check that gene intersection is non-empty before calling CuBlock.');
    end
    if n_probes < k
        error('kmeans: n_probes (%d) < k (%d). Cannot cluster fewer probes than requested clusters.', n_probes, k);
    end
    ...
```

Fix C alone is not sufficient (the Python guard in Fix A is required), but it makes debugging much easier if a 0-gene matrix ever reaches Octave for other reasons.

---

## 4. Issue 3 — Reduced Gene Coverage for All P/Q Datasets (MODERATE)

### What was observed

| Job | Common genes | Input genes | Coverage |
|---|---|---|---|
| `ANTE_Q0std` | 2366 | 3447 | 68.6% |
| `ANTE_QNBKass` | 3447 | 3447 | 100.0% |
| `NBBags_Q0std` | 2366 | 3447 | 68.6% |
| `GTExAffy_*` | 0 | 3447 | 0% (wrong orientation) |
| P0_standard (May 15 log) | 2293 | 3447 | 66.5% |

### Root Cause

The S0 strict dataset has 3447 genes — the no-NA intersection across ALL 88 cohorts and all platforms (Affymetrix, Illumina microarray, Illumina NGS, Agilent). This is an extremely restrictive filter: genes that weren't measured on even one platform are excluded.

The Zenodo reference P/Q datasets (P0, ANTE, NBBags, etc.) were designed for the original Shambhala2 paper on a different gene universe, primarily focused on microarray platforms. The ~31% loss for Q0std-based variants means 1154 genes from the strict FL gene set are not in Q0std.

The `ANTE_QNBKass` pair achieves 100% coverage because `Normal_B_Kassandra.csv` (QNBKass) appears to be an FL-project-specific reference dataset whose gene set was built to cover the same 3447 genes.

**The gene loss is largest when `Q0std` is used as the Q reference.** `QNBKass` avoids the loss entirely for ANTE. This is the key insight for selecting best P/Q combinations.

### Consequence

With `strategy='drop'` NA handling, the 1081 genes that are not in the intersection are dropped before harmonization and restored as NaN in the output. The resulting TSV uploaded to S3 will have ~31% NaN columns. This could distort the PCA R² batch metrics computed downstream.

### Investigation actions

1. For each P dataset, log the intersection size and fraction when generating calibration datasets.
2. Check whether `QNBKass` achieves better coverage for all P variants (not just ANTE) — if so, `QNBKass` should be the preferred Q for the strict S0 analysis.
3. For knn/softimpute S0 variants: these datasets have many more genes (up to ~21k instead of 3447), so the gene intersection with P/Q will be much larger. Run those variants to compare coverage.

---

## 5. Issue 4 — knn NA Strategy Consideration

The user asks: **"Maybe there are hidden NAs in the input dataframe, could `--na-strategy knn` help?"**

### Analysis

For the **strict** imputation S0:
- `NA strategy 'drop': 0 genes dropped.` — confirmed in all 5 observed jobs.
- The strict S0 dataset was built specifically to have zero NAs across all 7174 samples × 3447 genes. No hidden NAs are present in the expression values.
- `--na-strategy knn` would have no effect: there are no NaN genes to impute.

For the **knn** and **softimpute** S0 variants:
- These datasets have a larger gene set (genes with some NAs, now imputed). Shambhala's `apply_na_strategy` is set to `drop` in `run_shambhala_job.py` with the explicit rationale: *S0 knn/softimpute datasets are already imputed — double-imputing distorts values*.
- `--na-strategy knn` should NOT be used for these variants.

### Where hidden NAs could actually exist — Calibration P and Q datasets

Neither `run_shambhala_job.py` nor `run_shambhala.py` checks P or Q for NaN after downloading. If a calibration CSV has NaN values (e.g., from probe mapping gaps, empty cells, or `"NA"` strings not recognized as floats by pandas), those NaN values would:
1. Survive the gene intersection step (NaN is a value, not a missing column name)
2. Be piped to Octave as-is
3. Corrupt the quantilenorm → CuBlock computation silently (CuBlock states explicitly: "The micro array cannot contain NaN values")

**Recommended addition**: after downloading P and Q, assert they are NaN-free:

```python
if p_df.isna().any().any():
    n_na_cols = p_df.isna().any(axis=0).sum()
    raise ValueError(f"P calibration dataset '{P_DATASETS[p_key]}' has NaN values in {n_na_cols} gene columns. Fix the calibration CSV.")

if q_df.isna().any().any():
    n_na_cols = q_df.isna().any(axis=0).sum()
    raise ValueError(f"Q calibration dataset '{Q_DATASETS[q_key]}' has NaN values in {n_na_cols} gene columns.")
```

This does not affect the current crashes (the GTExAffy crash is caused by the empty intersection, not NaN), but it prevents a silent correctness bug for any P/Q dataset that has sparse NaN.

---

## 6. Issue 5 — DtypeWarning from Annotation Loading (MINOR)

### Symptom

```
DtypeWarning: Columns (5,6,8,...) have mixed types. Specify dtype option on import or set low_memory=False.
```

Emitted by `download_ann_from_s3` at `shambhala_bench_shared.py:121`. This is a warning, not an error.

### Cause

The annotation TSV has clinical columns (OS, PFS, cohort labels, etc.) that contain a mix of numeric values, strings, and empty cells. Pandas' default low-memory chunked parsing assigns different types to different chunks, then warns. This does not affect the expression data.

### Fix

Add `low_memory=False` to the annotation read in `download_ann_from_s3`:

```python
# Before
return pd.read_csv(gz, sep="\t", index_col=0)

# After
return pd.read_csv(gz, sep="\t", index_col=0, low_memory=False)
```

---

## 7. Implementation Plan

### Priority 1 — Fix GTExAffy calibration dataset (before next run)

**Where:** S3 bucket `your-bucket`, key `FL_batch_correction/calibration_datasets/GTExAffymetrix.csv`  
**Who does it:** Daniil must re-prepare and re-upload this file.  
**Steps:**
1. Download current file and inspect orientation.
2. Transpose to FL convention (samples × genes) and check that column names are HGNC symbols.
3. If columns are probe IDs, apply GPL570 → HGNC mapping (same logic as `GPL570_mapping.py`).
4. Save as `GTExAffymetrix_FL_convention.csv` locally, verify shape with `df.shape`, then overwrite on S3.
5. Run the smoke test to confirm intersection > 0: `python harmonization_scripts/run_shambhala_job.py --imp strict --method shambhala_GTExAffy_Q0std --n-shambhala-workers 1 --octave-bin /usr/bin/octave --random-seed 42 --out-json /tmp/gtex_test.json`

### Priority 2 — Add empty intersection guard (code fix, low risk)

**Files to edit:** `run_shambhala_job.py` and `run_shambhala.py`

In `run_shambhala_job.py`, after the existing log line at ~line 228:

```python
_log(label,
     f"Gene intersection: {len(common_genes)} common "
     f"(dropped {n_dropped_input} from input, {n_dropped_p} from P, {n_dropped_q} from Q)")

# ADD THIS BLOCK:
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
```

In `run_shambhala.py`, add after line 235 (after logging common gene count):

```python
if len(common_genes) == 0:
    logger.error(
        "Gene intersection is empty. Input, P, and Q must share gene names. "
        "Check that all CSVs use HGNC gene symbols in FL convention (rows=samples, cols=genes)."
    )
    sys.exit(1)
```

**Also add defensive guard to `kmeans.m`** (lines 5–8):

```matlab
[n_probes, n_samples] = size(data);
if n_probes == 0
    error('kmeans: data has 0 rows. Gene intersection was empty — check Python-level guard.');
end
if n_probes < k
    error('kmeans: n_probes=%d < k=%d. Reduce k or increase gene count.', n_probes, k);
end
```

### Priority 3 — Add NaN assertions for calibration P and Q

**File:** `run_shambhala_job.py`  
**Where:** After downloading P and Q (after line ~212 where `p_df` and `q_df` are assigned)

```python
# After downloading p_df and q_df:
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
```

Same pattern for `run_shambhala.py` after reading P and Q (around line 213–219).

### Priority 4 — Fix DtypeWarning in annotation loading

**File:** `harmonization_scripts/shambhala_bench_shared.py`  
**Change:** Line 121, in `download_ann_from_s3`:

```python
# Before
return pd.read_csv(gz, sep="\t", index_col=0)
# After
return pd.read_csv(gz, sep="\t", index_col=0, low_memory=False)
```

### Priority 5 — Log gene intersection as fraction (observability)

**File:** `run_shambhala_job.py`  
**Change:** Improve the intersection log message to show fraction:

```python
coverage_pct = 100 * len(common_genes) / max(len(input_T.index), 1)
_log(label,
     f"Gene intersection: {len(common_genes)} common ({coverage_pct:.1f}% of input) "
     f"(dropped {n_dropped_input} from input, {n_dropped_p} from P, {n_dropped_q} from Q)")
```

---

## 8. Expected Outcome After Fixes

| Scenario | Before fix | After fix |
|---|---|---|
| GTExAffy variants (6 jobs) | Octave crash, exit=1 | Either clear Python error if CSV still wrong, or runs successfully if CSV fixed |
| Any variant with 0-gene intersection | Octave crash with cryptic `r_sort` error | Clear Python error with diagnostic message |
| Any calibration CSV with hidden NaN | Silent NaN propagation in Octave output | Clear ValueError before dispatching to Octave |
| DtypeWarning | Noisy log output per job | Suppressed |
| ANTE / NBBags / P0 variants | Run with 66–69% gene coverage | No change (gene loss is a data issue, not a code bug) |
| GTExAffy variants after CSV fix | (as above) | Run with expected gene coverage |

---

## 9. Notes on Current Run State

The `shambhala_run_260516.log` was captured early (197 lines, ~5 minutes of a multi-hour run). The dispatcher is not crashed. The remaining 49 jobs (minus 2 already failed GTExAffy variants) are continuing to run.

To check current S3 output count from the Mac:
```bash
aws s3 ls s3://your-bucket/FL_batch_correction/exp/ | grep shambhala | wc -l
```

The GTExAffy failures are logged in `failed_jobs_shambhala.txt` inside the pod at `/app/Shambhala_containerized/harmonization_scripts/`. After fixing the CSV, re-run just the failed jobs:
```bash
python harmonization_scripts/run_shambhala_parallel.py \
    --n-workers 3 --n-shambhala-workers 5 \
    --octave-bin /usr/bin/octave --skip-if-exists --retry-failed \
    --memory-limit-gb 10.0 --timeout-s 21600 --random-seed 42
```

---

## 10. To-Do List — Code Changes

Tasks are ordered by priority. Each item maps to the implementation section above.

### Guard: empty gene intersection

- [x] **`run_shambhala_job.py`** — after the intersection log line (~line 228), add `if len(common_genes) == 0:` block that prints a diagnostic message to stderr, writes a sidecar with `status="empty_intersection"`, and calls `sys.exit(1)`. *(Priority 2)*
- [x] **`run_shambhala.py`** — after line 235 (after logging common gene count), add the same `if len(common_genes) == 0:` guard that logs the error and calls `sys.exit(1)`. *(Priority 2)*
- [x] **`octave/kmeans.m`** — add a size check at the top of the function body (after `size(data)`): error if `n_probes == 0`, error if `n_probes < k`. *(Priority 2, defensive)*

### Guard: NaN in calibration P and Q

- [x] **`run_shambhala_job.py`** — after downloading `p_df` and `q_df` (line ~212), add `p_df.isna().any().any()` and `q_df.isna().any().any()` checks that raise `ValueError` with the filename and column count. *(Priority 3)*
- [x] **`run_shambhala.py`** — after reading P and Q (lines ~213–219), add the same NaN assertions for `p_df` and `q_df`. *(Priority 3)*

### Logging: gene coverage fraction

- [x] **`run_shambhala_job.py`** — replace the plain intersection log line with one that also shows `coverage_pct = 100 * len(common_genes) / max(len(input_T.index), 1)` formatted as `({coverage_pct:.1f}% of input)`. *(Priority 5)*

### Suppress DtypeWarning

- [x] **`harmonization_scripts/shambhala_bench_shared.py`** — in `download_ann_from_s3`, add `low_memory=False` to `pd.read_csv`. *(Priority 4)*

---

## 11. Summary Table

| Issue | Severity | Root Cause | Fix |
|---|---|---|---|
| GTExAffy P dataset: 0 gene intersection | Critical | CSV stored in genes×samples (not FL convention) | Transpose + re-upload to S3 |
| No empty-intersection guard | Critical | Missing Python check before Octave dispatch | Add guard in `run_shambhala_job.py` and `run_shambhala.py` |
| kmeans crashes on 0-gene input | Critical (consequence of above) | Octave `rand(0,1)` → empty index | Add defensive check in `kmeans.m` |
| P/Q not checked for NaN | Moderate | No assertion after `download_calib_df` | Add NaN assertions after P/Q download |
| 31–34% gene loss for Zenodo P/Q refs | Moderate (data quality) | S0 strict gene set ≠ Zenodo gene universe | Use QNBKass or knn/softimpute S0 (more genes) |
| DtypeWarning from annotation loading | Minor | `low_memory=False` missing in pandas read | Add `low_memory=False` |
| `--na-strategy knn` | Not applicable | Strict S0 has 0 NaN; P/Q NaN is separate issue | Not needed for observed crashes; add P/Q NaN check instead |
