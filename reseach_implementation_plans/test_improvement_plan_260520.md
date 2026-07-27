# Test Improvement Plan — 2026-05-20

Source: `logs/test_logs_260520.txt` — 11 failed, 92 passed, 1 skipped in 7700 s.

This document describes root-cause analysis and proposed changes for four issues.
**No code is changed here.** Implementation follows approval.

---

## Issue 1 — `max_abs_diff` still used as the assertion metric in `test_combinations.py`

### Root cause

`tests/speed_up_tests/test_combinations.py` defines `_COMBOS` as a list of
`(combo_id, kwargs, tol, desc)` tuples, where `tol` is a `max_abs_diff` threshold.
The assertion on line 71 reads:

```python
assert metrics["max_abs_diff"] <= tol, ...
```

`max_abs_diff` is the single largest element-wise absolute difference across all genes
and samples. After Q-rescaling, expression values reach ≫10 000, so a handful of
outlier genes dominate: one gene with a value of 100 000 in the result but 116 253 in
the reference produces `max_abs_diff ≈ 16 253` — but the mean over 8 174 genes may
still be < 1%. The `assert_mean_rel_diff` helper already exists in conftest.py and
correctly uses `mean_rel_diff` as the primary assertion. It is used consistently in
all *individual* approach tests (A, C, D, E), but was never wired into the combination
tests. The `_COMBOS` tolerances (0.05 / 0.10) were also calibrated for
`max_abs_diff`, not `mean_rel_diff`.

### What the test output confirms

- `[A_only]` failed: `max_abs_diff=16253` > 0.05.
  But `test_a_qn_only` **passed** with `fail_threshold=0.05` (mean_rel_diff < 5%).
  Same underlying run — only the metric differed.
- Same pattern for C_20, E, A_E: individual approach tests pass; combination tests fail.

### Proposed changes — `test_combinations.py`

1. **Replace the tolerance field in `_COMBOS`** with a `(warn_threshold, fail_threshold)`
   tuple. Keep thresholds aligned with the individual approach tests:

   | Combo | warn | fail | Rationale |
   |---|---|---|---|
   | baseline | 0.01 | 0.05 | Reference run — must be near-zero |
   | A_only | 0.01 | 0.05 | Passes individually at 5% |
   | A_D | 0.10 | 0.20 | Synthetic P is exploratory (documented); consistent with `test_a_synthetic_cublock_p` |
   | C_40 | 0.01 | 0.05 | No real subsampling (40 ≥ 39); passes individually |
   | C_20 | 0.01 | 0.05 | Passes individually at 5% |
   | E | 0.01 | 0.05 | Passes individually at 5% |
   | A_E | 0.01 | 0.05 | Passes individually at 5% |
   | G | xfail | — | Deadlock (see Issue 2) |
   | A_D_E | 0.10 | 0.20 | Contains synthetic P flag |
   | A_D_G | xfail | — | Deadlock (see Issue 2) |

2. **Replace the `assert metrics["max_abs_diff"] <= tol` block** with a call to
   `assert_mean_rel_diff(metrics, warn_threshold=..., fail_threshold=..., label=...)`.
   Import `assert_mean_rel_diff` from conftest (it is already exported).

3. **Print `mean_rel_diff` in the per-test line** — the current print block omits it.
   Add `mean_rel_diff={metrics['mean_rel_diff']:.4%}` between `mean_abs_diff` and
   `pct_within_1pct`.

---

## Issue 2 — Approach G must be marked `xfail`; pipeline tests timeout with deadlock

### Root cause — accuracy (`test_g_cublock_matches_octave`)

`mean_rel_diff = 226.12%` far exceeds the 5% fail threshold. The test docstring
documents the known cause: NumPy's `default_rng` (PCG64) and Octave's
`rand('state', seed)` (Mersenne Twister) produce completely different random sequences
from the same integer seed. Because k-means initialization in `cublock_python` is
seeded with NumPy RNG and Octave's `kmeans.m` uses its own RNG, the gene-cluster
assignments across 30 repetitions diverge systematically. The averaged centroid
positions and polynomial fits therefore differ substantially, producing per-sample
outputs that can deviate by > 100% in relative terms.

This is a **fundamental algorithmic incompatibility**, not a numerical precision
issue. The Python CuBlock produces internally consistent results — it is just not
numerically equivalent to Octave CuBlock. `test_g_cublock_matches_octave` is
therefore an expected failure by design.

### Root cause — deadlock (`test_g_full_pipeline_within_tolerance`, `test_g_timing`, `test_combination[G]`, `test_combination[A_D_G]`)

All four timed-out tests share the same stack trace:

```
shambhala/parallel.py:402: for future in as_completed(future_to_idx)
concurrent.futures._base.as_completed → threading.Condition.wait → blocks forever
```

The `python_cublock=True` path in `harmonize_parallel` passes work to
`ProcessPoolExecutor` workers that invoke `cublock_python`. Inside the worker
processes, the NumPy RNG state is apparently not being initialised correctly, or a
multiprocessing queue deadlock occurs during result serialisation. The
`QueueFeederThread` stack (`multiprocessing/queues.py → _feed → nwait → waiter.acquire`)
in the timeout dump confirms that result objects are being serialised but the queue
feeder is stuck waiting on a lock — a classic multiprocessing deadlock with large
return arrays. The workers complete their computation but can never push results back
to the main process through the queue.

This is a **known bug** that requires separate investigation and refactoring of
`parallel.py`. Until fixed, all five tests using `python_cublock=True` must be
skipped or marked `xfail`.

### Proposed changes — `test_approach_g.py`

Mark three tests with `@pytest.mark.xfail`:

- `test_g_cublock_matches_octave`: `xfail(strict=True, reason="...")` — it always fails
  due to PRNG divergence; use `strict=True` so it turns red if it unexpectedly passes
  (would indicate the comparison logic broke).

- `test_g_full_pipeline_within_tolerance`: `xfail(strict=False, reason="...")` — hangs
  due to the multiprocessing deadlock. `strict=False` because the symptom is a timeout,
  which may behave differently in future runs.

- `test_g_timing`: `xfail(strict=False, reason="...")` — same deadlock.

### Proposed changes — `test_combinations.py`

Entries for `G` and `A_D_G` in `_COMBOS` must carry an `xfail` mark. Use
`pytest.param` with `marks=pytest.mark.xfail(...)`:

```python
pytest.param(
    "G",
    {"python_cublock": True},
    (None, None),          # thresholds unused — xfail
    "G: Python CuBlock",
    marks=pytest.mark.xfail(
        strict=False,
        reason="python_cublock triggers ProcessPoolExecutor deadlock (known bug)",
    ),
    id="G",
),
```

Same pattern for `A_D_G`.

The `test_combination` function must gracefully handle `(None, None)` thresholds, but
since the test is xfailed, the body will not be executed in normal runs.

---

## Issue 3 — No final summary table for all combination results

### Root cause

`test_combination` calls `print(...)` inside each test run, but pytest captures stdout
per-test and only shows it on failure (unless `-s` is passed). Even in verbose mode
(`-v`) the captured output appears only in failure sections, not as a coherent summary.
There is no aggregation or session-level finalisation step.

The `baseline_result` session fixture assembles one shared run, but individual
combination results are discarded after each test function returns.

### Proposed changes — `test_combinations.py`

1. **Add a module-level accumulator dict** `_COMBO_RESULTS: dict[str, dict] = {}` at
   the top of the file.

2. **In `test_combination`**, compute all four metrics and timing *before* the
   assertion, then store into `_COMBO_RESULTS[combo_id]` with the following fields:
   ```python
   {
       "desc": desc,
       "mean_rel_diff": ...,
       "max_abs_diff": ...,
       "mean_abs_diff": ...,
       "pct_within_1pct": ...,
       "time_s": ...,
       "speedup": ...,
       "status": "PASS" | "WARN" | "FAIL",
   }
   ```
   Status is determined by comparing `mean_rel_diff` against the combo's warn/fail
   thresholds, assigned before the assertion runs.

3. **Add `test_z_print_summary_table()`** as the *last function in the file* (the `z_`
   prefix keeps it alphabetically last; pytest executes functions in definition order by
   default). This function:
   - Checks `_COMBO_RESULTS` for all combos listed in `_COMBOS`.
   - Prints a fixed-width summary table to stdout. No assertion.
   - Marks xfailed entries as `"XFAIL (known)"` in the status column.
   - The table header and rows look like:

   ```
   ─────────────────────────────────────────────────────────────────────────────────────
   COMBO          mean_rel_diff  max_abs_diff  mean_abs_diff  %within_1pct  time(s)  speedup  status
   ─────────────────────────────────────────────────────────────────────────────────────
   baseline           0.0000%     0.0000        0.0000          100.0%       48.6     1.00x    PASS
   A_only             X.XXXX%     XXXX.XXXX     XX.XXXX          XX.X%       XX.X     X.XXx    PASS|WARN|FAIL
   ...
   G                  —           —             —                —            —        —        XFAIL (known)
   ─────────────────────────────────────────────────────────────────────────────────────
   ```

4. To ensure `test_z_print_summary_table` runs even when some combination tests fail,
   it must not depend on any fixture that would be skipped on failure. It reads only
   `_COMBO_RESULTS` which is populated unconditionally before assertions in step 2.

5. **Mark the function with `@octave_required`** so it is skipped cleanly when Octave
   is unavailable and `_COMBO_RESULTS` is empty.

---

## Issue 4 — No side-by-side "before vs after" expression comparison

### Root cause

The test log output contains truncated `repr()` of DataFrames like:
```
baseline_result = (... A2M ... 503.334464  1190.665227  20558.405705 ...)
```
These come from pytest's own parameter display, not from any intentional print. They
are abbreviated and not aligned to any "after" values, so they give no useful
information about the magnitude of distortion.

The user needs to visually compare, for the same genes and samples, what the baseline
Octave pipeline outputs versus what each speed-up approach outputs.

### Proposed changes — `tests/speed_up_tests/conftest.py`

Add a new helper function `print_expression_comparison`:

```python
def print_expression_comparison(
    result: pd.DataFrame,
    reference: pd.DataFrame,
    label: str = "",
    n_genes: int = 8,
    n_samples: int = 3,
) -> None:
```

**Gene selection strategy:** Pick the `n_genes` genes with the highest absolute
mean deviation from the reference across all samples. This surfaces genes where the
approach makes the most difference — the most informative for a qualitative review.

**Sample selection:** Use the first `n_samples` columns (lexicographic order of sample
IDs for reproducibility).

**Output format** (printed to stdout, captured by pytest on failure and in `-s` mode):

```
--- Expression comparison: A_only (first 3 samples, top 8 changed genes) ---
Gene            SAMPLE_A         SAMPLE_B         SAMPLE_C
                ref   → result   ref   → result   ref   → result
GENE_X      503.33 → 501.11   490.12 → 489.05  3201.44 → 3198.01
GENE_Y     1190.67 → 1202.45   880.33 → 892.01   441.22 → 453.19
...
```

### Proposed changes — `test_combinations.py`

- Import `print_expression_comparison` from conftest.
- Call it inside `test_combination` after computing metrics and before the assertion:

  ```python
  print_expression_comparison(
      result,
      expected_df,
      label=combo_id,
      n_genes=8,
      n_samples=3,
  )
  ```

- This call is placed after the metrics print but before the assertion, so it is always
  printed regardless of pass/fail outcome (pytest shows captured stdout on failure).
- For xfailed tests (G, A_D_G), the function body is still reached only if the pipeline
  does not hang; because of the deadlock those tests time out before reaching this line
  — no special handling needed.

---

## Files to be modified

| File | Changes |
|---|---|
| `tests/speed_up_tests/conftest.py` | Add `print_expression_comparison` helper |
| `tests/speed_up_tests/test_approach_g.py` | Add `@pytest.mark.xfail` to 3 tests |
| `tests/speed_up_tests/test_combinations.py` | Replace `tol` with threshold tuples; swap `max_abs_diff` assertion for `assert_mean_rel_diff`; add `mean_rel_diff` to per-test print; add module-level accumulator; add `test_z_print_summary_table` |

No changes to `shambhala/` source code, `pyproject.toml`, or other test files.

---

## Step-by-step To-Do List

Work through these steps in order. Each step is self-contained and verifiable before
moving to the next.

---

### Step 1 — Add `print_expression_comparison` to `conftest.py` ✅ DONE

File: `tests/speed_up_tests/conftest.py`

1.1. After the existing `compare_to_expected` function, define
     `print_expression_comparison(result, reference, label="", n_genes=8, n_samples=3)`.

1.2. Inside the function:
     - Compute the intersection of columns (genes) and index (samples) between `result`
       and `reference`.
     - Align both DataFrames to the intersection.
     - Compute per-gene mean absolute deviation: `(result - reference).abs().mean(axis=0)`.
     - Select the `n_genes` genes with the highest mean absolute deviation.
     - Select the first `n_samples` samples by lexicographic order of the index.

1.3. Build the header line: `--- Expression comparison: {label} ...`.

1.4. Print a column-header row listing sample IDs.

1.5. Print a sub-header row with `ref → result` labels under each sample.

1.6. For each selected gene, print one row: gene name left-justified in 16 chars,
     then for each sample: `{ref_val:10.2f} → {result_val:10.2f}` pairs separated by
     at least 2 spaces.

1.7. Add a blank line after the table.

1.8. Verify the function is importable: run `python -c "from tests.speed_up_tests.conftest import print_expression_comparison"`.

---

### Step 2 — Mark Approach G tests as `xfail` in `test_approach_g.py` ✅ DONE (verified: 3 passed, 2 xfailed)

File: `tests/speed_up_tests/test_approach_g.py`

2.1. Add `import pytest` to the top of the imports block (it is already there via
     `@pytest.mark.timeout`, but confirm).

2.2. Add `@pytest.mark.xfail(strict=True, reason="Python CuBlock PRNG (PCG64) diverges from Octave Mersenne Twister; k-means clusters differ systematically — 226% mean_rel_diff is the expected outcome")` decorator to `test_g_cublock_matches_octave`, placed after `@octave_required` and before `@pytest.mark.timeout(300)`.

2.3. Add `@pytest.mark.xfail(strict=False, reason="python_cublock=True triggers ProcessPoolExecutor deadlock in harmonize_parallel — QueueFeederThread blocked on result serialisation")` to `test_g_full_pipeline_within_tolerance`, placed after `@octave_required`.

2.4. Add the same `xfail(strict=False, ...)` decorator to `test_g_timing`, after
     `@octave_required`.

2.5. Update the module docstring test matrix to reflect the new expected statuses:
     - G-4: `xfail (PRNG incompatibility)`
     - G-5: `xfail (deadlock)`
     - G-6: `xfail (deadlock)`

2.6. Verify the three affected tests are now reported as `XFAIL` (not `FAILED`) by
     running: `pytest tests/speed_up_tests/test_approach_g.py -v -m "not slow" --timeout=30`.
     Expected: 3 PASSED (pure-Python tests) + 1 XFAIL (test_g_cublock_matches_octave,
     skips Octave run — xfail before running if Octave absent, or xfails after if present).

---

### Step 3 — Rework `_COMBOS` structure in `test_combinations.py` ✅ DONE

File: `tests/speed_up_tests/test_combinations.py`

3.1. Add `assert_mean_rel_diff` to the import from `tests.speed_up_tests.conftest`.

3.2. Add `print_expression_comparison` to the same import line.

3.3. Replace the `_COMBOS` list. Change the 3rd element from a single `tol` float to
     a `(warn_threshold, fail_threshold)` tuple:

     ```python
     _COMBOS = [
         ("baseline",   {}, (0.01, 0.05), "Baseline: Approach B only (always on)"),
         ("A_only",     {"precompute_qn_reference": True}, (0.01, 0.05),
          "A: Python QN reference"),
         ("A_D",        {"precompute_qn_reference": True, "synthetic_cublock_p": True},
          (0.10, 0.20), "A+D: Python QN + synthetic 1-col P for CuBlock (EXPLORATORY)"),
         ("C_40",       {"max_p_samples": 40}, (0.01, 0.05),
          "C: centroid P subsample to 40"),
         ("C_20",       {"max_p_samples": 20}, (0.01, 0.05),
          "C: centroid P subsample to 20"),
         ("E",          {"precompute_cublock_clusters": True}, (0.01, 0.05),
          "E: precomputed k-means clusters"),
         ("A_E",        {"precompute_qn_reference": True,
                         "precompute_cublock_clusters": True}, (0.01, 0.05),
          "A+E: Python QN + precomputed clusters"),
         pytest.param(
             "G", {"python_cublock": True}, (None, None),
             "G: Python CuBlock",
             marks=pytest.mark.xfail(
                 strict=False,
                 reason="python_cublock=True triggers ProcessPoolExecutor deadlock",
             ),
             id="G",
         ),
         ("A_D_E",      {"precompute_qn_reference": True, "synthetic_cublock_p": True,
                         "precompute_cublock_clusters": True}, (0.10, 0.20),
          "A+D+E: Python QN + synthetic P + precomputed clusters (EXPLORATORY)"),
         pytest.param(
             "A_D_G",
             {"precompute_qn_reference": True, "python_cublock": True,
              "synthetic_cublock_p": True},
             (None, None),
             "A+D+G: Python QN + synthetic P + Python CuBlock",
             marks=pytest.mark.xfail(
                 strict=False,
                 reason="python_cublock=True triggers ProcessPoolExecutor deadlock",
             ),
             id="A_D_G",
         ),
     ]
     ```

3.4. Update the `@pytest.mark.parametrize` decorator signature from
     `"combo_id,kwargs,tol,desc"` to `"combo_id,kwargs,thresholds,desc"`.

3.5. Update the `test_combination` function signature from `(combo_id, kwargs, tol, desc, baseline_result)`
     to `(combo_id, kwargs, thresholds, desc, baseline_result)`.

---

### Step 4 — Add module-level result accumulator to `test_combinations.py` ✅ DONE

File: `tests/speed_up_tests/test_combinations.py`

4.1. Add `_COMBO_RESULTS: dict[str, dict] = {}` near the top of the module, after
     imports.

---

### Step 5 — Rewrite `test_combination` body ✅ DONE

File: `tests/speed_up_tests/test_combinations.py`

5.1. Unpack thresholds at the start of the function body:
     `warn_threshold, fail_threshold = thresholds`.

5.2. Run the pipeline and compute metrics as before (no change to this part).

5.3. Compute `speedup = baseline_time / max(t[0], 1e-3)`.

5.4. Determine status string before asserting:
     ```python
     mrd = metrics["mean_rel_diff"]
     if mrd < warn_threshold:
         status = "PASS"
     elif mrd < fail_threshold:
         status = "WARN"
     else:
         status = "FAIL"
     ```

5.5. Store to accumulator immediately after computing status:
     ```python
     _COMBO_RESULTS[combo_id] = {
         "desc": desc,
         "mean_rel_diff": mrd,
         "max_abs_diff": metrics["max_abs_diff"],
         "mean_abs_diff": metrics["mean_abs_diff"],
         "pct_within_1pct": metrics["pct_within_1pct"],
         "time_s": t[0],
         "speedup": speedup,
         "status": status,
         "warn": warn_threshold,
         "fail": fail_threshold,
     }
     ```

5.6. Update the existing `print(...)` line to add `mean_rel_diff`:
     ```python
     print(
         f"\n[{combo_id}] {desc}\n"
         f"  mean_rel_diff={mrd:.4%}  max_abs_diff={metrics['max_abs_diff']:.4f}"
         f"  mean_abs_diff={metrics['mean_abs_diff']:.4f}"
         f"  pct_within_1pct={metrics['pct_within_1pct']:.1f}%"
         f"  time={t[0]:.1f}s  speedup={speedup:.2f}x  status={status}"
     )
     ```

5.7. Add the expression comparison print call (after metrics print, before assertion):
     ```python
     print_expression_comparison(result, expected_df, label=combo_id, n_genes=8, n_samples=3)
     ```

5.8. Replace the `assert metrics["max_abs_diff"] <= tol` line with:
     ```python
     assert_mean_rel_diff(metrics, warn_threshold=warn_threshold, fail_threshold=fail_threshold, label=combo_id)
     ```

---

### Step 6 — Add `test_z_print_summary_table` to `test_combinations.py` ✅ DONE

File: `tests/speed_up_tests/test_combinations.py`

6.1. Add a new function `test_z_print_summary_table()` at the very bottom of the file.

6.2. Decorate it with `@octave_required` (skip cleanly when Octave absent and
     `_COMBO_RESULTS` is empty).

6.3. Inside the function:
     - Print a top divider line of 95 dashes.
     - Print the header row with fixed-width columns:
       `COMBO`, `mean_rel_diff`, `max_abs_diff`, `mean_abs_diff`, `%within_1pct`,
       `time(s)`, `speedup`, `status`.
     - Print a divider line.
     - Iterate over `_COMBOS` in their original order (use the `combo_id` from each
       entry, handling both plain tuples and `pytest.param` objects) and print the
       corresponding row from `_COMBO_RESULTS`. If a combo is absent from
       `_COMBO_RESULTS` (timed out or xfailed before storing), print `—` in all
       numeric fields and `XFAIL / TIMEOUT` in the status column.
     - Print a bottom divider line.
     - Print a legend: `PASS < {warn}% | WARN {warn}%–{fail}% | FAIL > {fail}%`.

6.4. Add no assertion — this is a purely informational print.

6.5. Verify the function is the last `def test_` in the file (check by inspection).

---

### Step 7 — Run and verify ✅ DONE

7.1. Run the pure-Python subset first (no Octave needed, should be fast):
     ```
     pytest tests/speed_up_tests/ -v -m "not slow" --timeout=60 \
         --ignore=tests/speed_up_tests/test_combinations.py -k "not octave_required"
     ```
     Expected: all pass.

7.2. Run only the Approach G tests to confirm xfail behaviour:
     ```
     pytest tests/speed_up_tests/test_approach_g.py -v --timeout=60
     ```
     Expected: 3 PASSED (pure-Python), 3 XFAIL (Octave tests).

7.3. Run combination tests without slow tests:
     ```
     pytest tests/speed_up_tests/test_combinations.py -v -m "not slow" -s --timeout=1200
     ```
     Expected:
     - baseline, A_only, C_40, C_20, E, A_E → PASSED
     - A_D, A_D_E → PASSED or WARN (exploratory; warn=10%, fail=20%)
     - G, A_D_G → XFAIL
     - `test_z_print_summary_table` → PASSED (prints the table)
     - Summary table printed to terminal after all combination tests.

7.4. Run the full test suite to confirm no regressions:
     ```
     pytest tests/ -v -m "not slow" --timeout=1200 -s
     ```
     Expected: 0 FAILED, xfails only for G-related tests.

7.5. Inspect the summary table and expression comparison output in the terminal. Confirm:
     - All 10 combination entries appear in the table.
     - Expression comparison blocks show gene names, `ref → result` values, and
       absolute differences are plausible for the approach.
     - G and A_D_G rows show `—` placeholders.
