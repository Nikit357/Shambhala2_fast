# Test Failures Analysis — 2026-05-14

**Test run:** `pytest tests/ -v` from `Shambhala_containerized/`
**Results:** 8 failed, 22 passed
**Environment:** Octave is installed at `/opt/conda/bin/octave`

---

## Summary of Failures

| Test | Error class | Root cause |
|---|---|---|
| `test_bridge_returns_correct_shape` | `OctaveExecutionError` | Bug 1 — `source()` relative path |
| `test_bridge_symbol_index_correct` | `OctaveExecutionError` | Bug 1 — `source()` relative path |
| `test_single_worker_matches_reference` | `OctaveExecutionError` → `RuntimeError` | Bug 1 — `source()` relative path |
| `test_parallel_matches_single` | `OctaveExecutionError` → `RuntimeError` | Bug 1 — `source()` relative path |
| `test_local_io_roundtrip` | `OctaveExecutionError` → `RuntimeError` | Bug 1 — `source()` relative path |
| `test_gene_intersection_is_applied` | `OctaveExecutionError` → `RuntimeError` | Bug 1 — `source()` relative path |
| `test_na_knn_full_pipeline` | `OctaveExecutionError` → `RuntimeError` | Bug 1 + Bug 2 (samples dropped before Octave) |
| `test_na_drop_full_pipeline` | `ValueError: number sections must be larger than 0` | Bug 2 — samples dropped instead of genes |

After fixing Bug 1, Bugs 2 and 3 will cause the integration tests to fail with wrong results or empty output. All three bugs must be fixed.

---

## Bug 1 — `source('Shambhala2_piped.m')` uses a relative path

### Location
`shambhala/octave_bridge.py`, line 203 (inside `run_octave_normalize`):
```python
eval_preamble = (
    f"addpath('{octave_scripts_dir}'); "
    f"{seed_preamble}"
    f"NH={nh}; NP={np_}; k={k}; "
    f"source('Shambhala2_piped.m');"   # ← BUG
)
```

### Root cause
`addpath(octave_scripts_dir)` adds the `octave/` directory to Octave's **function search path**, which allows calls like `CuBlock(...)` and `quantilenorm(...)` to resolve. However, **`source()` resolves relative paths from Octave's current working directory** (inherited from the Python subprocess, which is wherever pytest or the CLI was invoked from — e.g. the `Shambhala_containerized/` root). The file `Shambhala2_piped.m` is in `octave/`, not in the cwd, so Octave fails with:

```
error: no such file, '/...Shambhala_containerized/Shambhala2_piped.m'
```

### Evidence
Octave stderr in every failing test:
```
error: no such file, '/home/jovyan/.../Shambhala_containerized/Shambhala2_piped.m'
```
The path shown is the pytest working directory, not `octave/`.

### Fix
Replace the relative `source()` with the absolute path:
```python
f"source('{octave_scripts_dir}/Shambhala2_piped.m');"
```
One-line change. No other files need to change for this bug.

### Tests affected
All 8 failing tests. `test_na_drop_full_pipeline` is an exception — it crashes before reaching Octave due to Bug 2, but would also hit this bug if Bug 2 were not present.

---

## Bug 2 — `na_handling.apply_na_strategy` and `restore_na_genes` called with genes×samples instead of samples×genes

### Location
**`run_shambhala.py`**, Step 5 (line ~247) and Step 8 (line ~286):
```python
# Step 5 (wrong orientation)
clean_T, na_mask = na_handling.apply_na_strategy(input_T, ...)  # input_T is genes×samples

# Step 8 (wrong orientation)
full_harmonized_T = na_handling.restore_na_genes(harmonized_T, na_mask)
# harmonized_T is genes×samples; restore_na_genes expects samples×genes
```

**`tests/test_integration.py`**, `_run_pipeline` helper, lines ~54 and ~68:
```python
clean_T, na_mask = na_handling.apply_na_strategy(input_T, strategy=na_strategy)
...
full_harmonized_T = na_handling.restore_na_genes(harmonized_T, na_mask)
```

### Root cause
The `na_handling` module was designed for **FL convention (samples×genes)**. Its docstring specifies `df : pd.DataFrame, shape (n_samples, n_genes)`. All logic inside assumes:
- **Rows = samples**, **Columns = genes**

In `detect_na_genes`:
```python
return list(df.columns[df.isna().any(axis=0)])  # columns = genes in FL convention
```

But both callers transpose before calling:
```python
input_T = input_df.T  # now rows=genes, columns=samples
clean_T, na_mask = na_handling.apply_na_strategy(input_T, ...)
```

When `input_T` is passed (genes×samples):
- `df.columns` = **sample names** (not gene names)
- `detect_na_genes(input_T)` returns the **sample IDs** that have any NaN gene, not the genes
- `df.drop(columns=na_genes)` drops those **sample columns** instead of gene rows
- For `input_with_nas.csv` (where every sample has at least one injected NaN), all 10 samples are detected and dropped
- `clean_T.shape = (8174 genes, 0 samples)` → 0 columns

In `harmonize_parallel`:
```python
n_samples = input_df.shape[1]    # = 0
effective_workers = min(1, 0)    # = 0
np.array_split([], 0)            # ValueError: number sections must be larger than 0
```

The `NAMask.original_gene_order` also stores **sample names** instead of gene names, so `restore_na_genes` would produce wrong output even if it didn't crash.

The **strategy='knn'** path has the same bug: `na_frac = df.isna().mean(axis=0)` computes fraction per column (per sample), so it drops samples (not genes) that exceed `max_na_frac`.

### Evidence
```
# test_na_drop_full_pipeline traceback:
shambhala/parallel.py:164: in harmonize_parallel
    for arr in np.array_split(sample_names, effective_workers)
ary = [], indices_or_sections = 0, axis = 0
ValueError: number sections must be larger than 0.
```

### Fix
In both `run_shambhala.py` and `tests/test_integration.py`:

1. Pass `input_T.T` (back to samples×genes) to `apply_na_strategy`; then re-transpose the cleaned result:
   ```python
   clean_df, na_mask = na_handling.apply_na_strategy(input_T.T, strategy=na_strategy, ...)
   clean_T = clean_df.T  # genes×samples, NaN-free, for Octave
   ```

2. Filter `p_T` to only the genes that survived NA drop (otherwise `pd.concat` in `_harmonize_one_batch` would fill missing genes with NaN, corrupting Octave input):
   ```python
   p_T = p_T.loc[clean_T.index]
   ```

3. Pass `harmonized_T.T` to `restore_na_genes` and stop transposing the result in Step 9 (it's already samples×genes):
   ```python
   result_df = na_handling.restore_na_genes(harmonized_T.T, na_mask)  # samples×genes
   # Step 9 no longer needs .T — result_df is already FL convention
   ```

### Tests affected
`test_na_drop_full_pipeline` crashes today (0 samples). `test_na_knn_full_pipeline` crashes due to Bug 1 first, but would fail via the same mechanism if Bug 1 were fixed.

---

## Bug 3 — `compute_q_statistics` called with genes×samples instead of samples×genes

### Location
**`run_shambhala.py`**, Step 6 (line ~261):
```python
rm, rs = q_rescale.compute_q_statistics(q_T, q_pseudocount=args.q_pseudocount)
# q_T = q_df.T — genes×samples, but compute_q_statistics expects samples×genes
```

**`tests/test_integration.py`**, `_run_pipeline`, line ~55:
```python
rm, rs = q_rescale.compute_q_statistics(q_T)
```

### Root cause
`compute_q_statistics` is implemented for **FL convention (samples×genes)**:
```python
gene_symbols = list(q_df.columns)                          # columns = genes in FL convention
rm = pd.Series(np.mean(log_q, axis=0), index=gene_symbols) # axis=0 = per-column = per-gene
```

When called with `q_T` (genes×samples):
- `q_T.columns` = **GTEX sample IDs** (not gene names)
- `gene_symbols` = 100 GTEX sample names
- `np.mean(log_q, axis=0)` computes mean across 8174 gene rows for each of 100 sample columns → 100 per-sample means
- `rm` and `rs` are indexed by GTEX sample names, length 100

### Evidence
Printed in every integration test traceback:
```
rm = GTEX.14JIY.1426.SM.6EU29    7.168095
     GTEX.YF7O.2426.SM.5IFJL     7.181746
     ...
     Length: 100, dtype: float64
```
Should be 8174 gene-symbol values; instead 100 GTEX sample names.

### Downstream effect
`rescale(octave_output, rm, rs)` computes:
```python
common_genes = octave_output.index.intersection(rm.index)
```
`octave_output.index` = 8174 gene names. `rm.index` = 100 GTEX sample names. Intersection = **empty set**. `rescale` returns an empty DataFrame (0 rows). `harmonize_parallel` then returns empty output for all tests.

This bug is **masked by Bug 1** today (tests fail at Octave before reaching `rescale`), but will cause all integration tests to produce empty output once Bug 1 is fixed.

### Fix
In both `run_shambhala.py` and `tests/test_integration.py`, pass the genes filtered to survivors of NA drop, transposed back to samples×genes:
```python
q_df_filtered = q_T.loc[clean_T.index].T  # samples×genes, filtered to NA-surviving genes
rm, rs = q_rescale.compute_q_statistics(q_df_filtered, q_pseudocount=args.q_pseudocount)
```

The filtering step (`q_T.loc[clean_T.index]`) also aligns Q to exactly the genes that Octave will output, avoiding downstream gene-set mismatches.

### Tests affected
All 6 integration tests (would fail with empty output after Bug 1 is fixed).

---

## Implementation Plan

All three bugs must be fixed together. Bugs 2 and 3 are in the same two files at adjacent lines, so they are fixed in the same pass.

### Task 1 — Fix `source()` path in `octave_bridge.py`

**File:** `shambhala/octave_bridge.py`

Change line 203:
```python
# Before
f"source('Shambhala2_piped.m');"

# After
f"source('{octave_scripts_dir}/Shambhala2_piped.m');"
```

**Validation:** `test_bridge_returns_correct_shape` and `test_bridge_symbol_index_correct` should pass.

---

### Task 2 — Fix orientation mismatch in `run_shambhala.py`

**File:** `run_shambhala.py`

Steps 5, 6, 8, and 9 of `main()` need to change. The gene intersection (Step 4) and the transpose step (Step 3) stay the same.

**Step 5 — Apply NA strategy (current / fixed):**
```python
# Before
clean_T, na_mask = na_handling.apply_na_strategy(input_T, ...)

# After
clean_df, na_mask = na_handling.apply_na_strategy(input_T.T, ...)
clean_T = clean_df.T
p_T = p_T.loc[clean_T.index]  # align P to surviving genes
```

**Step 6 — Compute Q statistics (current / fixed):**
```python
# Before
rm, rs = q_rescale.compute_q_statistics(q_T, q_pseudocount=args.q_pseudocount)

# After
q_df_filtered = q_T.loc[clean_T.index].T  # samples×genes, filtered to surviving genes
rm, rs = q_rescale.compute_q_statistics(q_df_filtered, q_pseudocount=args.q_pseudocount)
```

**Step 8 — Restore NA genes (current / fixed):**
```python
# Before
full_harmonized_T = na_handling.restore_na_genes(harmonized_T, na_mask)
# ...
# Step 9
result_df = full_harmonized_T.T

# After
result_df = na_handling.restore_na_genes(harmonized_T.T, na_mask)
# result_df is already samples×genes; Step 9 (.T) is removed
```

---

### Task 3 — Fix orientation mismatch in `tests/test_integration.py`

**File:** `tests/test_integration.py`

The `_run_pipeline` helper has the same orientation bugs as `run_shambhala.py`.

**Lines ~54–69 (current):**
```python
clean_T, na_mask = na_handling.apply_na_strategy(input_T, strategy=na_strategy)
rm, rs = q_rescale.compute_q_statistics(q_T)
harmonized_T = parallel_module.harmonize_parallel(input_df=clean_T, p_df=p_T, ...)
full_harmonized_T = na_handling.restore_na_genes(harmonized_T, na_mask)
return full_harmonized_T.T
```

**Fixed:**
```python
clean_df, na_mask = na_handling.apply_na_strategy(input_T.T, strategy=na_strategy)
clean_T = clean_df.T
p_T = p_T.loc[clean_T.index]
q_df_filtered = q_T.loc[clean_T.index].T
rm, rs = q_rescale.compute_q_statistics(q_df_filtered)
harmonized_T = parallel_module.harmonize_parallel(input_df=clean_T, p_df=p_T, ...)
result_df = na_handling.restore_na_genes(harmonized_T.T, na_mask)
return result_df
```

---

### To-do list

- [x] **Task 1**: Fix `source()` path in `shambhala/octave_bridge.py` (1 line)
- [x] **Task 2**: Fix orientation of `apply_na_strategy`, `compute_q_statistics`, `restore_na_genes` calls in `run_shambhala.py` (Steps 5, 6, 8, 9)
- [x] **Task 3**: Fix orientation of the same calls in `tests/test_integration.py` (`_run_pipeline`)
- [x] **Task 4**: Run full test suite and confirm 30/30 pass (or 22 pass + 8 skip if Octave is unavailable — but since Octave IS installed here, all 30 should pass)
- [x] **Task 5**: (Optional) Add `p_T` gene-alignment to the `harmonize_parallel` function signature documentation to make the gene-alignment contract explicit

### Files to change

| File | Lines changed | Bug fixed |
|---|---|---|
| `shambhala/octave_bridge.py` | 1 line (source path) | Bug 1 |
| `run_shambhala.py` | ~6 lines (Steps 5, 6, 8, 9) | Bugs 2 + 3 |
| `tests/test_integration.py` | ~6 lines (`_run_pipeline`) | Bugs 2 + 3 |

No changes needed to `na_handling.py`, `q_rescale.py`, `parallel.py`, or any `.m` files.
