# Plan: Add `I_rare_batches_removed` Strategy to Shambhala Cross-Product Benchmark

**Date:** 2026-05-19  
**Author:** Daniil Nikitin  
**Status:** ✅ COMPLETE — implemented 2026-05-19

---

## Background

The main harmonization benchmark (`harmonization-scripts/bench_shared.py`) received a new filter strategy `I_rare_batches_removed` (commit `cde2a87`, 2026-05-16). It removes all RNA batches with fewer than **50 samples**, keeping ~5,000–5,300 of the ~6,500 S0 samples across 15 removed batches (including `RNASeq_FFPE_PolyA`, `RNASeq_FF_Total`, and 13 microarray/NGS minority batches).

The Shambhala cross-product scripts (`harmonization_scripts/` inside `Shambhala_containerized/`) still list only 11 strategies in `shambhala_bench_shared.py:ALL_STRATEGIES`; strategy I is absent. This plan describes adding it.

---

## Key Property: Shambhala Is Sample-Independent

Shambhala normalizes each sample independently. The harmonized value for sample X in `S0_no_removal__{imp}__{method}__post0.tsv.gz` is **identical** to what you would get by running Shambhala on a dataset that contains only strategy I samples. The samples do not affect each other's normalization.

This means:

> If `S0_no_removal__{imp}__{method}__post0.tsv.gz` already exists on S3, the strategy I output can be derived by **filtering rows** — no re-harmonization required.

This avoids relaunching all 54 Shambhala jobs (which take hours each on K8s) just to add a new strategy.

---

## Two Complementary Paths

### Path 1 — Fast derivation from existing S0 outputs (primary, runs first) ✅ COMPLETE

A new script `harmonization_scripts/derive_rare_batches_outputs.py`:

1. Downloads `S0_no_removal__strict__ann.tsv.gz` from `prepared/` (uses strict imputation; strategy I sample IDs do not depend on imputation since the rare-batch filter operates purely on the `RNA_BATCH` column).
2. Computes strategy I sample IDs: keep samples from batches where `ann[RNA_BATCH].value_counts() >= 50`. This replicates the `MIN_BATCH_SIZE = 50` logic from `bench_shared.py:84`.
3. Uploads the filtered annotation under all three imputation keys so that `run_shambhala_job.py` can find them in future runs:
   - `prepared/I_rare_batches_removed__strict__ann.tsv.gz`
   - `prepared/I_rare_batches_removed__knn__ann.tsv.gz`
   - `prepared/I_rare_batches_removed__softimpute__ann.tsv.gz`
4. For each `(imp, method)` pair (3 × 18 = 54 combinations):
   a. Check whether `S0_no_removal__{imp}__{method}__post0.tsv.gz` exists on S3.
   b. If not: skip this pair (will be handled by Path 2 when the full job runs).
   c. If yes:
      - Download `S0_no_removal__{imp}__{method}__post0.tsv.gz`.
      - Filter rows to strategy I sample IDs → `exp_I`.
      - Upload as `I_rare_batches_removed__{imp}__{method}__post0.tsv.gz`.
      - Apply `identify_outlier_batches(exp_I, ann_I)`, remove outlier batches → `exp_I_post1`.
      - Upload as `I_rare_batches_removed__{imp}__{method}__post1.tsv.gz`.
      - Record metrics row: `r2_batch`, `r2_diag`, `n_samples`, `n_genes`, `status="derived"`.
5. After all pairs: appends derived rows to the existing sidecar JSONs at `shambhala_metrics/{imp}__{method}.json` (or creates them if absent), then re-uploads `shambhala_metrics.csv`.

**Note on post1:** The `S0_no_removal__{imp}__{method}__post1.tsv.gz` file must **not** be used as input for strategy I's post0 — it already has S0-level outlier batches removed, which differs from I-level outliers. Always start from `post0` and apply `identify_outlier_batches` fresh on the strategy-I filtered subset.

**Command:**
```bash
python harmonization_scripts/derive_rare_batches_outputs.py \
    --skip-if-exists \
    --log-level INFO
```

**CLI arguments to support:**
- `--skip-if-exists`: check S3 before downloading / uploading (default should be True in practice)
- `--imps`: restrict to specific imputation(s) (default: all 3)
- `--methods`: restrict to specific method(s) (default: all 18)
- `--log-level`: DEBUG / INFO / WARNING

---

### Path 2 — Integration into existing jobs (for future re-runs) ✅ COMPLETE

After Path 1 runs, all outputs exist on S3. The existing `run_shambhala_job.py` and `run_shambhala_parallel.py` still do not know about strategy I. This matters for future re-runs (e.g., if a job fails and needs to be retried, or if new P/Q calibration datasets are added).

Changes required:

#### `harmonization_scripts/shambhala_bench_shared.py`

- **Line 51–63:** Add `"I_rare_batches_removed"` to `ALL_STRATEGIES`:
  ```python
  ALL_STRATEGIES: list[str] = [
      "S0_no_removal",
      "A_confirmed_bad",
      ...
      "H_affymetrix_extended",
      "I_rare_batches_removed",   # ← add
  ]
  ```
  This automatically propagates to every loop in `run_shambhala_job.py` (which iterates `ALL_STRATEGIES`) and to `_all_outputs_exist` in `run_shambhala_parallel.py` (which checks `ALL_STRATEGIES × 2` outputs per job).

#### `harmonization_scripts/run_shambhala_job.py`

No logic changes are required. The `for strat in ALL_STRATEGIES` loop at lines 381 and 389 already handles any strategy generically — it downloads the annotation from `prepared/{strat}__{imp}__ann.tsv.gz` (uploaded by Path 1) and filters the harmonized matrix. The only change is updating the docstring and comment counts:

- **Line 8:** `{strat}__{imp}__ann.tsv.gz   (×11 strategies)` → `(×12 strategies)`
- **Line 3:** `22 (strategy × post_rm) outputs` → `24 (strategy × post_rm) outputs`
- **Line 13:** `(×22)` → `(×24)`
- **Line 398:** log message `[{n_done}/22]` → `[{n_done}/24]`

#### `harmonization_scripts/run_shambhala_parallel.py`

- **Line 6:** `Each worker writes 22 (strategy × post_rm) outputs → 1,188 total S3 objects.` → `24 outputs → 1,296 total S3 objects.`
- **Line 45 (docstring):** `all 22 outputs` → `all 24 outputs`
- **Line 228 (docstring of `_all_outputs_exist`):** `all 22 output keys` → `all 24 output keys`

No logic changes — `_all_outputs_exist` already derives the expected keys from `ALL_STRATEGIES` dynamically.

---

## Execution Sequence

Run on the K8s pod `shambhala-bench` (16 vCPU / 128 GiB):

```bash
# Step 1: Run fast derivation (uses existing S0 outputs; ~minutes per method)
python harmonization_scripts/derive_rare_batches_outputs.py \
    --skip-if-exists \
    --log-level INFO \
    > /workspace/derive_rare_batches.log 2>&1

# Step 2 (optional, only if some (imp, method) pairs were skipped in Step 1):
# Run the full job dispatcher with --skip-if-exists so it only processes missing outputs.
# After shambhala_bench_shared.py is updated with strategy I, the dispatcher will
# launch 54 jobs. Each job will skip all 24 outputs that already exist and only
# produce the missing ones.
nohup python harmonization_scripts/run_shambhala_parallel.py \
    --n-workers 3 \
    --n-shambhala-workers 5 \
    --octave-bin /usr/bin/octave \
    --skip-if-exists \
    --memory-limit-gb 10.0 \
    --random-seed 42 \
    > /workspace/shambhala_fill_missing.log 2>&1 &
```

Step 2 is only needed if the previous full run did not complete all 54 jobs (some had `status=3` — S0 not found). In a typical scenario where all 1,188 existing outputs are present, Step 1 alone produces all 108 new strategy I outputs (54 × 2 post_rm).

---

## Output Count Delta

| Before | After | Delta |
|---|---|---|
| 11 strategies | 12 strategies | +1 |
| 22 outputs per job | 24 outputs per job | +2 |
| 1,188 total S3 exp objects | 1,296 total S3 exp objects | +108 |
| 54 sidecar JSONs (22 rows each) | 54 sidecar JSONs (24 rows each) | +108 rows |

The 108 new keys follow the existing pattern:
```
FL_batch_correction/exp/I_rare_batches_removed__{imp}__{method}__post{0|1}.tsv.gz
```
where `imp` ∈ {strict, knn, softimpute} and `method` ∈ 18 Shambhala variants.

---

## Files to Create or Modify ✅ ALL COMPLETE

| File | Action | Status |
|---|---|---|
| `harmonization_scripts/derive_rare_batches_outputs.py` | **Create** | ✅ Created |
| `harmonization_scripts/shambhala_bench_shared.py` | **Edit** | ✅ `"I_rare_batches_removed"` added to `ALL_STRATEGIES`; `[union-attr]` → `[attr-defined]` |
| `harmonization_scripts/run_shambhala_job.py` | **Edit** | ✅ Docstring counts updated (22→24, ×11→×12); `[union-attr]` → `[attr-defined]` |
| `harmonization_scripts/run_shambhala_parallel.py` | **Edit** | ✅ Docstring counts updated (22→24, 1188→1296); `[union-attr]` → `[attr-defined]` |
| `harmonization_scripts/k8s/README.md` | **Edit** | ✅ 1,188 → 1,296 in pod description and monitoring target |
| `harmonization_scripts/README.md` | **Edit** | ✅ 11→12 strategies, 1,188→1,296 in table, monitoring target, and S3 paths section; derive script quickstart added |
| `CLAUDE.md` (this directory) | **Edit** | ✅ 1,188→1,296, 11→12 strategies, 22→24 outputs per job; derive script added to file table |

---

## Open Questions / Assumptions

1. **S0 annotation availability:** The derivation script assumes `prepared/S0_no_removal__strict__ann.tsv.gz` exists on S3. This is created by `run_prep_parallel.py` in `harmonization-scripts/` and should exist since the main benchmark ran. If it is missing for some imputation, the script should fall back to the next available imputation (all give the same sample IDs for strategy I).

2. **Prepared annotation for strategy I:** The main benchmark's `run_prep_parallel.py` now includes `I_rare_batches_removed` (added in commit `cde2a87`), but it may not have been run yet. The derivation script creates these annotation files as a side effect (Step 3 of Path 1), so `run_shambhala_job.py` will find them in future runs.

3. **`MIN_BATCH_SIZE` coupling:** The derivation script hardcodes `MIN_BATCH_SIZE = 50` inline (same constant as `bench_shared.py:84`). If the threshold changes in the future, both files must be updated together. Consider extracting this constant to `shambhala_bench_shared.py` if a more formal coupling is desired.

4. **Sidecar JSON deduplication:** The derivation script appends rows to existing sidecar JSONs. The `_upload_consolidated_metrics` function in `run_shambhala_parallel.py` deduplicates by `(strat, imp, method, post_rm)`, keeping the last entry. Running the derivation script and then a full dispatcher run will produce duplicates only for strategy I rows — this is handled correctly by the existing deduplication logic.
