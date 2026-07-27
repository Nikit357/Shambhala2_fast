# S0 Harmonized Output Reuse — Implementation Plan

**Date:** 2026-05-21  
**Author:** Daniil Nikitin  
**Status:** Draft — awaiting approval before implementation

---

## Problem Statement

`run_shambhala_job.py` currently has two modes:

1. **Full run** (default): download S0 prepared data + calibration P/Q, run Octave harmonization (~3–15 h per job depending on variant), write 24 outputs (12 strategies × 2 post_rm) to S3.
2. **Full cache** (`--skip-if-exists`, all 24 keys present): skip everything immediately.

There is no middle path. When new batch-removal strategies are added (e.g., strategy `I_rare_batches_removed` was added after the initial 54-job run), every job must re-run the full Octave harmonization even though the underlying harmonized matrix is identical. This wastes hours of compute.

**Key insight:** The S0 harmonized file — `FL_batch_correction/exp/S0_no_removal__{imp}__{method}__post0.tsv.gz` — already contains every harmonized sample. All other 23 outputs are row-subsets (sample-filtered views) of this matrix, optionally with outlier batches removed. If S0 is on S3, Octave never needs to run again for that `(imp × method)` pair.

This feature generalizes the pattern already used by `derive_rare_batches_outputs.py`.

---

## Three-Case Decision Tree

After Step 0 S3 checks, the job falls into one of three cases:

```
skip_if_exists=True?
├── NO  → Case C: full Octave run (existing behavior, no change)
└── YES → check all 24 output keys on S3
           ├── ALL present → Case A: early exit "cached" (existing behavior, no change)
           ├── SOME missing, S0_no_removal__post0 PRESENT → Case B: reuse path (NEW)
           └── S0_no_removal__post0 MISSING → Case C: full Octave run
```

| Case | Condition | Action | Octave invoked? |
|------|-----------|--------|-----------------|
| A | `--skip-if-exists` + all 24 keys exist | Early exit | No |
| B | `--skip-if-exists` + S0 harmonized exists + ≥1 missing | Download S0 file, derive missing | **No** |
| C | S0 harmonized missing OR `--skip-if-exists` absent | Full harmonization | **Yes** |

---

## Changes to `run_shambhala_job.py`

### 1. Refactor Step 0 (Early-Exit and Reuse Check)

**Current Step 0** (lines 193–203): checks if all 24 keys exist → exits. Only activates with `--skip-if-exists`.

**New Step 0**: same guard on `--skip-if-exists`, but instead of a single all-or-nothing check, it:

1. Builds the full list of 24 `(strat, post_rm, s3_key)` tuples.
2. Checks each key — records which are missing.
3. If none missing → Case A: early exit.
4. If some missing, checks `s3_key_exp("S0_no_removal", imp, method, False)` (the anchor key).
5. If anchor exists → Case B: calls `_run_reuse_path(...)`.
6. Otherwise → falls through to Case C (full harmonization, unchanged).

```python
# ── Step 0: Cache and reuse check ─────────────────────────────────────────────
if args.skip_if_exists:
    all_output_specs = [
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
        # Case A
        _log(label, "All 24 outputs already on S3 — cached.")
        rows = _make_failed_rows(imp, method, status="cached")
        _write_sidecar(rows, s3, imp, method, args.out_json)
        sys.exit(0)

    s0_anchor_key = s3_key_exp("S0_no_removal", imp, method, False)
    if s3_exists(s3, s0_anchor_key):
        # Case B
        _log(
            label,
            f"S0 harmonized found on S3 ({s0_anchor_key}). "
            f"Deriving {len(missing_specs)} missing output(s) without re-harmonizing.",
        )
        rows = _run_reuse_path(
            s3=s3,
            imp=imp,
            method=method,
            missing_specs=missing_specs,
            all_output_specs=all_output_specs,
            args=args,
            label=label,
        )
        _write_sidecar(rows, s3, imp, method, args.out_json)
        sys.exit(0)

    # Case C: S0 not on S3, fall through to full harmonization
    _log(label, f"S0 harmonized not found on S3 — running full harmonization.")
```

---

### 2. New `_run_reuse_path` Function

Add a new top-level function in `run_shambhala_job.py`. It mirrors Steps 5–8 of the existing `main()` but without any Octave involvement.

```python
def _run_reuse_path(
    s3: object,
    imp: str,
    method: str,
    missing_specs: list[tuple[str, bool, str]],
    all_output_specs: list[tuple[str, bool, str]],
    args: argparse.Namespace,
    label: str,
) -> list[dict]:
    """
    Derive missing (strat × post_rm) outputs from the pre-existing S0 harmonized
    matrix on S3. No Octave run — pure download, filter, upload.

    Parameters
    ----------
    missing_specs
        List of (strat, post_rm, s3_key) tuples for outputs not yet on S3.
    all_output_specs
        Complete list of all 24 (strat, post_rm, s3_key) tuples.

    Returns
    -------
    list[dict]
        One row per (strat × post_rm) pair for the sidecar JSON.
    """
```

**Internal steps:**

#### Step B1 — Memory pre-flight

Same check as the existing Step 1 — verify `free_gb() >= args.memory_limit_gb` before downloading the (potentially large) S0 harmonized matrix.

#### Step B2 — Download S0 harmonized matrix

```python
s0_key = s3_key_exp("S0_no_removal", imp, method, False)
_log(label, f"Downloading S0 harmonized matrix from s3://{S3_BUCKET}/{s0_key} ...")
t0 = time.time()
exp_harmonized = download_exp_from_s3(s3, s0_key)
_log(label, f"S0 harmonized: {exp_harmonized.shape[0]} samples × {exp_harmonized.shape[1]} genes ({time.time()-t0:.0f}s)")
```

This is identical to `download_exp_from_s3` used for prepared data — it returns a samples × genes DataFrame.

#### Step B3 — Download strategy annotations (missing strategies only)

Only download annotations for strategies that appear in `missing_specs`. Do not download all 12 unconditionally.

```python
strats_needed = {strat for strat, _, _ in missing_specs}
ann_strat: dict[str, pd.DataFrame] = {}
for strat in strats_needed:
    ann_strat[strat] = download_ann_from_s3(s3, s3_key_prepared_ann(strat, imp))
_log(label, f"Downloaded annotations for {len(strats_needed)} strategy/ies: {sorted(strats_needed)}")
```

#### Step B4 — Derive and upload each missing output

Loop only over `missing_specs`. For each:

```python
rows: list[dict] = []
n_done = 0
for strat, post_rm, key in missing_specs:
    n_done += 1
    tag = f"{strat}×post_rm={post_rm}"

    # Filter harmonized matrix to this strategy's sample IDs
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
        _log(label, f"  Upload complete.")
    except Exception:
        status = "upload_failed"
        _log(label, f"  Upload FAILED:\n{traceback.format_exc()}")

    rows.append({
        "strat": strat, "imp": imp, "method": method, "post_rm": post_rm,
        "r2_batch": r2_b, "r2_diag": r2_d,
        "n_samples": len(ann_s), "n_genes": exp_strat.shape[1],
        "status": status,
    })
    del exp_strat, ann_s
    gc.collect()
```

#### Step B5 — Append cached rows for already-existing outputs

After writing missing outputs, build placeholder rows for the outputs that already existed on S3, so the sidecar JSON is complete (all 24 entries).

```python
missing_keys = {key for _, _, key in missing_specs}
for strat, post_rm, key in all_output_specs:
    if key not in missing_keys:
        rows.append({
            "strat": strat, "imp": imp, "method": method, "post_rm": post_rm,
            "r2_batch": float("nan"), "r2_diag": float("nan"),
            "n_samples": float("nan"), "n_genes": float("nan"),
            "status": "cached",
        })
```

Return `rows` (24 entries total, sorted is fine).

---

## Changes to `run_shambhala_parallel.py`

**No changes required.** The dispatcher is transparent to this feature — it launches `run_shambhala_job.py` subprocesses and reads their exit codes and JSON sidecars. Case B exits with code 0, identical to Case A.

The only behavioral difference visible to the dispatcher: Case B jobs will take seconds-to-minutes instead of hours, and the sidecar will contain a mix of `"ok"` and `"cached"` statuses. Both are already handled correctly in `_run_dispatcher`.

Optional enhancement to the dispatcher (low priority, not required for correctness): in `_scan_existing_outputs`, also detect which jobs are Case B candidates (S0 present, some missing) and log them separately. This helps in progress monitoring. Can be added later.

---

## Changes to `shambhala_bench_shared.py`

**No new functions required.** The anchor key is `s3_key_exp("S0_no_removal", imp, method, False)`, which already exists. The `download_exp_from_s3` helper is reused as-is.

---

## Affected CLI Behavior

### `run_shambhala_job.py`

| Flag combination | Old behavior | New behavior |
|-----------------|--------------|--------------|
| No `--skip-if-exists` | Full run always | Full run always (unchanged) |
| `--skip-if-exists`, all 24 cached | Early exit | Early exit (unchanged) |
| `--skip-if-exists`, S0 present, ≥1 missing | N/A (was full run) | **Reuse path (new)** |
| `--skip-if-exists`, S0 missing, ≥1 missing | N/A (was full run) | Full run (unchanged) |

### `run_shambhala_parallel.py`

No behavioral change. Pass `--skip-if-exists` as before to activate reuse.

---

## Edge Cases and Design Decisions

### Gene set consistency

The S0 harmonized file was written with a specific gene set: the intersection of (input genes ∩ P genes ∩ Q genes) after NA drop. Any new strategy outputs derived from it will have exactly the same gene columns. This is correct — strategy filtering is row-only (sample IDs), never column-level.

### S0_no_removal post1 (outlier-batch-removed S0)

`S0_no_removal__post1` is a valid output — it applies `identify_outlier_batches` to the full S0 dataset. The reuse path derives it correctly by filtering rows of the downloaded `exp_harmonized` (the post0 S0 file). There is no circular dependency.

### Speed-up flags (precompute_qn_reference, etc.)

Case B never calls Octave, so speed-up flags are irrelevant. The downloaded harmonized matrix is the final result regardless of how it was originally produced. No special handling needed.

### Double-downloading the S0 anchor key

In Case B, we download `S0_no_removal__post0` as `exp_harmonized`. If `S0_no_removal__post0` itself is in `missing_specs` (i.e., it was missing when the check happened but somehow appeared on S3 just before download — race condition), the upload in Step B4 will overwrite it with identical content. This is benign and has no correctness impact.

### Partial S0 failure (S0 anchor present but corrupted)

If `download_exp_from_s3` raises an exception, wrap Step B2 in a try/except that falls back to Case C (full harmonization) with a clear log message. This prevents a corrupted S0 anchor from permanently blocking reuse.

### `--skip-if-exists` absent + S0 present

No reuse check is performed. The job runs full harmonization and overwrites all 24 outputs including S0. This is intentional — the absence of `--skip-if-exists` means the user wants a clean rerun.

### Memory footprint of Case B vs Case C

Case B downloads: S0 harmonized (~7200 × 3500 float64 ≈ 200 MB) + strategy annotations (negligible).  
Case C downloads: S0 prepared (same size) + P/Q calibration (~40 MB) + runs Octave (peaks at ~2× matrix size).  
Case B peak RAM is lower than Case C. The existing `memory_limit_gb` pre-flight is still appropriate.

---

## Testing Plan

### Unit tests (new file: `tests/test_s0_reuse.py`)

1. **`test_reuse_path_derives_correct_outputs`**  
   Mock S3 with a pre-populated `S0_no_removal__post0` file (real fixture data). Set missing_specs to include one non-S0 strategy (e.g. `A_confirmed_bad`). Call `_run_reuse_path`. Assert the uploaded DataFrame is the correct row-subset of the S0 harmonized matrix.

2. **`test_reuse_path_post_rm_applied`**  
   Same setup but `post_rm=True`. Patch `identify_outlier_batches` to return one batch name. Assert uploaded DataFrame has fewer rows.

3. **`test_step0_routes_case_a`**  
   All 24 S3 keys exist. Assert job exits with code 0 and sidecar has all `"cached"` statuses.

4. **`test_step0_routes_case_b`**  
   S0 anchor exists, two keys missing. Assert `_run_reuse_path` is called (mock it). Assert Octave is never invoked.

5. **`test_step0_routes_case_c_no_skip_flag`**  
   `--skip-if-exists` absent, S0 present. Assert full harmonization path is taken (Octave invoked).

6. **`test_step0_routes_case_c_no_s0`**  
   `--skip-if-exists` set, S0 anchor missing. Assert full harmonization path is taken.

### Integration test (manual, on pod)

1. Run one `(imp × method)` job to completion (all 24 outputs on S3).
2. Delete 3 specific strategy outputs.
3. Re-run same job with `--skip-if-exists`.
4. Confirm: only 3 outputs uploaded; Octave not invoked (no Octave log lines); sidecar has 3 `"ok"` + 21 `"cached"` entries.
5. Verify r2_batch of the re-derived outputs matches original values (within float precision).

---

## Implementation Checklist

- [x] Refactor Step 0 in `run_shambhala_job.py` into three-case logic (Case A / B / C).
- [x] Implement `_run_reuse_path` function in `run_shambhala_job.py`.
- [x] Add fallback try/except around S0 download in `_run_reuse_path` → fall through to Case C on failure.
- [x] Update docstring at top of `run_shambhala_job.py` to describe all three cases.
- [x] Add unit tests in `tests/test_s0_reuse.py` (15 tests, all pass).
- [ ] Run manual integration test on pod.
- [x] Update `harmonization_scripts/CLAUDE.md` — add Case B description to the exit-codes table and the architecture section.

---

## Files Modified

| File | Change type |
|------|------------|
| `harmonization_scripts/run_shambhala_job.py` | Refactor Step 0; add `_run_reuse_path` |
| `tests/test_s0_reuse.py` | New test file |
| `harmonization_scripts/CLAUDE.md` | Documentation update |
| `harmonization_scripts/run_shambhala_parallel.py` | None |
| `harmonization_scripts/shambhala_bench_shared.py` | None |
