# Shambhala Octave Timeout Failure Analysis — 2026-05-17

**Run command:**
```
python harmonization_scripts/run_shambhala_parallel.py \
  --n-workers 5 --n-shambhala-workers 5 \
  --octave-bin /usr/bin/octave --skip-if-exists \
  --memory-limit-gb 10.0 --timeout-s 216000 --random-seed 42
```

**Observed outcome:** 80 out of 80 Octave batch calls timed out at exactly 7200 seconds. Zero jobs completed successfully.

---

## 1. Root Cause #1 — Software Bug: Timeout Not Propagated to Job Workers

### What happened

`run_shambhala_parallel.py` operates on a two-level timeout hierarchy:

| Level | Variable | Controls |
|---|---|---|
| **Outer** | `--timeout-s` → `proc.wait(timeout=timeout_s)` | How long the dispatcher waits for the entire `run_shambhala_job.py` subprocess |
| **Inner** | `--shambhala-timeout-s` → `subprocess.run(..., timeout=timeout_s)` in `octave_bridge.py` | How long each individual Octave batch call is allowed to run |

The `launch()` function in `run_shambhala_parallel.py` builds the `cmd` list (lines 256–269) that it passes to `subprocess.Popen`. It correctly forwards `--n-shambhala-workers`, `--octave-bin`, `--memory-limit-gb`, `--random-seed`, and `--skip-if-exists` to `run_shambhala_job.py`. **It does not pass `--shambhala-timeout-s`.**

As a result, the inner Octave timeout in every job subprocess defaults to 7200 seconds — the `argparse` default in `run_shambhala_job.py` line 102:

```python
parser.add_argument("--shambhala-timeout-s", type=int, default=7200, dest="shambhala_timeout_s")
```

This is then passed to `harmonize_parallel()` (line 282), which passes it to each `_harmonize_one_batch()` (parallel.py line 79), which passes it to `run_octave_normalize()` (octave_bridge.py line 212) as the `timeout` argument of `subprocess.run`.

### Evidence

- The dispatcher logs show: `Timeout: 216000s` — the outer timeout was set correctly to 60 hours.
- All 80 `subprocess.TimeoutExpired` lines in the log show: `timed out after 7200 seconds` — exactly the hardcoded inner default.
- Every failing Octave command has `NH=1435` (all five workers hit the same limit simultaneously).
- No job ran for more than ~7220 seconds (7200s Octave + ~20s overhead).

### The call chain

```
run_shambhala_parallel.py
  --timeout-s 216000          ← correctly sets proc.wait(timeout=216000)
  launch() builds cmd list
    [sys.executable, run_shambhala_job.py,
     --imp, --method, --out-json, --memory-limit-gb,
     --n-shambhala-workers, --octave-bin, --random-seed]
     # ← --shambhala-timeout-s is ABSENT

run_shambhala_job.py
  args.shambhala_timeout_s = 7200   ← argparse default, never overridden
  harmonize_parallel(..., timeout_s=7200, ...)

parallel.py → _harmonize_one_batch(..., timeout_s=7200)
octave_bridge.py → subprocess.run([octave ...], timeout=7200)
  # ← Octave killed after 7200s regardless of --timeout-s 216000
```

---

## 2. Root Cause #2 — Computational Overload: NH=1435 per Octave Batch

Even if the timeout were wired correctly, the per-batch computation time would need careful sizing. The current batching is far too coarse for interpreted Octave.

### Numbers

- S0 strict dataset: **7174 samples × 3447 genes**
- `--n-shambhala-workers 5`: samples split into 5 batches → **NH = 7174 / 5 = 1435 per batch**
- All P calibration sets range from NP=141 (NBKass) to NP=878 (NBext)
- Each Octave call processes a pool matrix of **~2366–3447 genes × (1435 + NP) ≈ 1576–2313 columns**

### What CuBlock does with NH=1435

CuBlock (the inner Octave algorithm) performs the following **30 times per call**, then averages:

1. **k-means clustering of genes** (k=5): input matrix is n_genes × (NH + NP). For NH=1435 and NP=733, this is 2366 genes × 2168 columns. The custom `kmeans.m` uses pure-random centroid initialization and iterates until convergence. In interpreted Octave, each k-means iteration scans all n_genes × n_columns values — O(k × n_genes × (NH + NP)) per iteration, repeated until convergence.
2. **Cubic polynomial fitting per cluster**: fits a degree-3 polynomial over NH + NP values per gene per cluster.
3. QN precedes CuBlock: quantile normalization of the full pool (n_genes × (NH + NP)).

All operations scale **linearly with NH + NP** for the polynomial fitting and **quadratically with NH + NP** for k-means distance computations (distance from each gene-vector to each centroid). Interpreted Octave runs these as nested M-file loops with no JIT compilation, making the constant factor large.

### Why 7200 seconds is not enough

For NH=1435 and typical NP=200–800:
- Pool matrix: ~2400 genes × 1600–2300 samples
- CuBlock: 30 × k_iter × k × 2400 × 2300 floating-point operations in Octave interpreted loops
- Observed: every batch exceeded 2 hours without completing

For comparison, the toy test with 10 samples completes in ~2 minutes. Scaling from NH=10 to NH=1435 is a 143× increase in sample count, making the pool ~10× larger in total matrix size (since NP dominates for small NH), and Octave's interpreted loops do not scale sub-linearly.

### The CPU contention problem

With `--n-workers 5` and `--n-shambhala-workers 5`, the pod runs up to **25 simultaneous Octave processes** on a 16-vCPU pod. Each Octave process is single-threaded, so they compete for CPU. This adds further latency on top of the already slow computation.

---

## 3. Secondary Issue: GTExAffy Gene Intersection = 0

The two `GTExAffy` P variants failed immediately (in 36 seconds) with `Gene intersection: 0 common`. These were the calibration CSVs stored in genes×samples orientation (the transposed FL convention bug). Daniil reported these as already corrected; they are excluded from further analysis here.

---

## 4. How Many Jobs Were Affected

| Category | Count | Observation |
|---|---|---|
| Timed-out Octave calls (7200s) | 80 | All NH=1435; all strict/knn/softimpute NB* variants |
| Failed immediately (gene intersection=0) | 2 | GTExAffy, corrected separately |
| Completed successfully | 0 | No job reached S3 upload phase |
| Run aborted mid-softimpute | — | Log ends during softimpute runs; not all 54 jobs finished |

The dispatcher was still processing `softimpute` jobs at the last log entry. All `strict` (18 variants × 5 batches = 90 Octave calls; 80 timeout errors plus 10 from GTExAffy early-exit) and all `knn` jobs failed. The `softimpute` run was ongoing.

---

## 5. Implementation Plan

### Fix 1: Wire `--timeout-s` to `--shambhala-timeout-s` (critical, ~5 lines)

**File:** `harmonization_scripts/run_shambhala_parallel.py`
**Location:** `launch()` function, `cmd` list construction (~line 256–269)

Add `--shambhala-timeout-s` to the cmd list so the inner Octave timeout tracks the outer job timeout. The inner timeout should be equal to or slightly less than the outer timeout (the batches all start simultaneously and run in parallel, so they finish within one timeout window):

```python
# After the existing --random-seed block:
cmd.extend(["--shambhala-timeout-s", str(timeout_s)])
```

This single change makes `--timeout-s` in the dispatcher correctly control both the outer job wait and the inner Octave subprocess limit.

Optionally, add a dedicated `--shambhala-timeout-s` CLI argument to the dispatcher for independent control:

```python
parser.add_argument(
    "--shambhala-timeout-s",
    type=int,
    default=None,
    dest="shambhala_timeout_s",
    help=(
        "Per-Octave-batch timeout in seconds. "
        "Defaults to --timeout-s if not set."
    ),
)
```

Then in `_run_dispatcher`, accept `shambhala_timeout_s: int | None` and pass `shambhala_timeout_s or timeout_s` to `cmd.extend(["--shambhala-timeout-s", str(...)])`.

### Fix 2: Reduce NH per Batch (required to complete within any reasonable timeout)

The core problem is NH=1435 per Octave call. Increasing `--n-shambhala-workers` reduces NH proportionally:

| `--n-shambhala-workers` | NH per batch (7174 samples) | Pool columns (NP≈400) | Estimated time |
|---|---|---|---|
| 5 (current) | 1435 | 1835 | >2 hours (timed out) |
| 16 | 449 | 849 | ~30–40 min (estimate) |
| 50 | 144 | 544 | ~8–12 min (estimate) |
| 100 | 72 | 472 | ~5–8 min (estimate) |

However, increasing n_shambhala_workers per job while keeping `--n-workers 5` multiplies the simultaneous Octave process count. The pod has 16 vCPU. Recommended configuration:

```
--n-workers 1 \
--n-shambhala-workers 16 \
--timeout-s 86400 \
--shambhala-timeout-s 86400
```

This gives:
- 16 Octave processes simultaneously (matches vCPU count)
- NH = 7174 / 16 ≈ 449 per batch
- Estimated ~35 min per Octave batch → under 1 hour per job
- Jobs run sequentially (n-workers=1), no CPU contention between jobs

If faster throughput is needed with sequential safety:
```
--n-workers 2 --n-shambhala-workers 8 --timeout-s 86400 --shambhala-timeout-s 86400
```
(16 total Octave processes, 2 jobs in parallel, each using 8 of 16 vCPU)

### Fix 3: Add Per-Batch Timing to Logs (diagnostic, low effort)

**File:** `shambhala/parallel.py`, `_harmonize_one_batch()`

After the `run_octave_normalize()` call, log the elapsed time and NH:
```python
import time
t0 = time.time()
result = run_octave_normalize(...)
elapsed = time.time() - t0
logger.info("Batch NH=%d NP=%d completed in %.0fs.", len(batch_sample_names), np_, elapsed)
return result
```

This makes future runs produce per-batch timing in the logs, enabling accurate timeout estimation before committing to a full 54-job run.

### Fix 4: Add a Smoke-Test Batch Before Full Run (optional, defensive)

Before dispatching all 54 jobs, run a single small-NH test to confirm Octave is fast enough for the chosen `--n-shambhala-workers`:

```python
# In run_shambhala_parallel.py main(), before _run_dispatcher:
if args.smoke_test:
    _run_smoke_test(n_shambhala_workers=args.n_shambhala_workers, ...)
```

Alternatively, this is already available via `pytest tests/` — but not with production-scale NH.

---

## 6. Recommended Re-Run Command

After applying Fix 1 and Fix 2:

```bash
python harmonization_scripts/run_shambhala_parallel.py \
    --n-workers 1 \
    --n-shambhala-workers 16 \
    --octave-bin /usr/bin/octave \
    --skip-if-exists \
    --memory-limit-gb 10.0 \
    --timeout-s 7200 \
    --random-seed 42 \
    > /workspace/shambhala_run_260517b.log 2>&1 &
```

With Fix 1 applied, `--timeout-s 7200` now flows to `--shambhala-timeout-s 7200` in the job. With NH≈449 per batch, each Octave call should complete well within 7200 seconds, leaving margin for safety. If the first job completes in under 30 minutes, scale up to `--n-workers 2`.

---

## 7. File Change Summary

| File | Change | Priority |
|---|---|---|
| `harmonization_scripts/run_shambhala_parallel.py` | Add `--shambhala-timeout-s` to `cmd` in `launch()` | **Critical — blocks any timeout control** |
| `harmonization_scripts/run_shambhala_parallel.py` | Add `--shambhala-timeout-s` CLI arg (optional, for independent control) | Medium |
| `shambhala/parallel.py` | Log per-batch NH, NP, elapsed time in `_harmonize_one_batch` | Low (diagnostic) |
| `harmonization_scripts/k8s/README.md` | Document recommended n_workers / n_shambhala_workers for 16-vCPU pod | Low (docs) |

---

## 8. Step-by-Step Implementation To-Do List

All action items below are approved. They are ordered by dependency: Steps 1–3 must be done before any re-run; Steps 4–5 are diagnostic improvements that can be applied in the same commit.

---

### Step 1 — Fix the timeout propagation bug in `run_shambhala_parallel.py` (CRITICAL) ✅ DONE

**File:** `harmonization_scripts/run_shambhala_parallel.py`

**Why first:** Without this fix, every subsequent run will continue silently capping each Octave batch at 7200 s regardless of `--timeout-s`, making `--timeout-s` effectively decorative at the job level.

**Sub-steps:**

1.1. Open `run_shambhala_parallel.py`. Locate the `_run_dispatcher()` function signature (line ~224). Add a new parameter `shambhala_timeout_s: int | None` between `timeout_s` and `memory_limit_gb`.

1.2. In `main()`, add a new argparse argument directly after the existing `--timeout-s` argument:
```python
parser.add_argument(
    "--shambhala-timeout-s",
    type=int,
    default=None,
    dest="shambhala_timeout_s",
    help=(
        "Per-Octave-batch timeout in seconds passed to run_shambhala_job.py "
        "--shambhala-timeout-s. Defaults to --timeout-s if not set."
    ),
)
```

1.3. In the `main()` call to `_run_dispatcher()`, pass the resolved value:
```python
shambhala_timeout_s=args.shambhala_timeout_s or args.timeout_s,
```

1.4. Inside `_run_dispatcher()`, thread `shambhala_timeout_s` through to the `launch()` closure. Update the closure signature (or capture via nonlocal) so `launch()` can read it.

1.5. In `launch()`, inside the `cmd` list construction, add after the `--random-seed` block:
```python
cmd.extend(["--shambhala-timeout-s", str(shambhala_timeout_s)])
```

1.6. Update the log line that prints job parameters (the `print(f"Shambhala jobs to run: ...")` in `_run_dispatcher`) to also show the inner Octave timeout:
```
f"... | Timeout (job): {timeout_s}s | Timeout (Octave batch): {shambhala_timeout_s}s | ..."
```

1.7. Update the module-level docstring at the top of the file: add `--shambhala-timeout-s N` to the documented argument list with a clear explanation that it is the per-Octave-batch limit (distinct from `--timeout-s` which is the per-job limit).

---

### Step 2 — Verify the wiring with a dry run (no Octave needed) ⏳ Run yourself in the pod

**Why:** Confirm that the new argument appears in the subprocess command before doing any pod work.

**Sub-steps:**

2.1. From the `Shambhala_containerized/` directory, run:
```bash
python harmonization_scripts/run_shambhala_parallel.py \
    --imps strict --methods shambhala_P0std_Q0std \
    --n-workers 1 --n-shambhala-workers 2 \
    --timeout-s 9999 --shambhala-timeout-s 8888 \
    --dry-run   # (if dry-run mode exists) OR inspect the printed cmd
```
If no dry-run mode exists, temporarily add `print(cmd)` before `subprocess.Popen` in `launch()`, run once, verify `--shambhala-timeout-s 8888` appears in the printed command, then remove the debug print.

2.2. Also verify the fallback: run without `--shambhala-timeout-s` and confirm the log shows `Timeout (Octave batch): 9999s` (inheriting from `--timeout-s`).

---

### Step 3 — Add per-batch timing instrumentation to `parallel.py` (diagnostic) ✅ DONE

**File:** `shambhala/parallel.py`, function `_harmonize_one_batch()`

**Why:** Without per-batch timing, the next run produces no evidence of whether NH=449 is fast enough. If the next run also times out, there will be no data to estimate the actual safe NH.

**Sub-steps:**

3.1. Add `import time` at the top of `parallel.py` (if not already present).

3.2. In `_harmonize_one_batch()`, wrap the `run_octave_normalize()` call with a timer:
```python
t0 = time.time()
result = run_octave_normalize(
    pool_df=pool_df,
    nh=len(batch_sample_names),
    np_=p_df.shape[1],
    k=k,
    octave_scripts_dir=octave_scripts_dir,
    octave_bin=octave_bin,
    timeout_s=timeout_s,
    random_seed=random_seed,
)
elapsed = time.time() - t0
logger.info(
    "Octave batch done: NH=%d, NP=%d, elapsed=%.0fs.",
    len(batch_sample_names), p_df.shape[1], elapsed,
)
return result
```

3.3. Confirm that `logger` is already set up in `parallel.py` (it is: `logger = logging.getLogger(__name__)` at line 19). No additional setup needed.

3.4. In `run_shambhala_job.py`, verify that the logging level is propagated so that `logger.info` messages from `parallel.py` appear in the job's stdout. If not, add before the `harmonize_parallel()` call:
```python
import logging
logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                    format="%(asctime)s %(name)s %(levelname)s %(message)s")
```

---

### Step 4 — Determine the correct `--n-shambhala-workers` value for the 16-vCPU pod ⏳ Run yourself in the pod

**Why:** The recommended value of 16 is an estimate. Before committing to a full 54-job run, validate with a single test job on the real pod.

**Sub-steps:**

4.1. rsync the updated code to the pod (standard procedure per `harmonization_scripts/k8s/README.md`).

4.2. Run a single test job with NH≈449 and measure elapsed time:
```bash
python harmonization_scripts/run_shambhala_job.py \
    --imp strict \
    --method shambhala_P0std_Q0std \
    --n-shambhala-workers 16 \
    --octave-bin /usr/bin/octave \
    --shambhala-timeout-s 7200 \
    --random-seed 42 \
    --out-json /tmp/test_job.json
```
Record the wall-clock time. Expected: 30–60 min for NH≈449.

4.3. If the test job finishes under 30 minutes, consider increasing to `--n-shambhala-workers 30` (NH≈239) for the full run to leave a comfortable time buffer.

4.4. If the test job finishes in 30–60 minutes, keep `--n-shambhala-workers 16` and set `--shambhala-timeout-s 7200` with headroom.

4.5. If the test job exceeds 60 minutes, increase to `--n-shambhala-workers 50` (NH≈144) and re-test.

4.6. Record the confirmed NH and elapsed time in the `harmonization_scripts/k8s/README.md` as a benchmark reference for future runs.

---

### Step 5 — Update `harmonization_scripts/k8s/README.md` with pod configuration guidance ✅ DONE

**File:** `harmonization_scripts/k8s/README.md`

**Why:** The current README does not document the relationship between `--n-workers`, `--n-shambhala-workers`, pod vCPU count, and NH. This led to the NH=1435 mistake. Future runs should have the guidance inline.

**Sub-steps:**

5.1. Add a new section "Choosing `--n-workers` and `--n-shambhala-workers`" explaining:
- Total simultaneous Octave processes = `n-workers × n-shambhala-workers`
- This should not exceed the pod's vCPU count (16 for the current shambhala pod)
- NH per batch = `n_S0_samples / n-shambhala-workers`; CuBlock scales roughly linearly with NH+NP
- Recommended for 16-vCPU pod: `--n-workers 1 --n-shambhala-workers 16` (NH≈449) as a safe starting point

5.2. Add a subsection "Observed timing data" with a table to be filled in after Step 4:
| n-shambhala-workers | NH (strict S0) | Elapsed per job | Date tested |
|---|---|---|---|
| 5 | 1435 | >7200s (timed out) | 2026-05-17 |
| 16 | 449 | TBD | — |

5.3. Add the two-level timeout hierarchy explanation (outer `--timeout-s` vs inner `--shambhala-timeout-s`) so it is clear to future operators which knob controls what.

---

### Step 6 — Full re-run on the pod ⏳ Run yourself in the pod

**Why:** Execute the corrected benchmark.

**Sub-steps:**

6.1. Confirm all previously completed S3 outputs are still intact (use `--skip-if-exists`). Since no jobs succeeded in the failed run, all 54 jobs need to be re-run. Remove or clear `harmonization_scripts/failed_jobs_shambhala.txt` if it was written during the failed run, so that `--skip-if-exists` (not `--retry-failed`) governs job selection.

6.2. Launch the corrected run using the timing-validated parameters from Step 4 (example for confirmed 16-workers case):
```bash
nohup python harmonization_scripts/run_shambhala_parallel.py \
    --n-workers 1 \
    --n-shambhala-workers 16 \
    --octave-bin /usr/bin/octave \
    --skip-if-exists \
    --memory-limit-gb 10.0 \
    --timeout-s 7200 \
    --random-seed 42 \
    > /workspace/shambhala_run_260517b.log 2>&1 &
```

6.3. Monitor the first job in the log. Confirm that:
- The per-batch timing lines appear (from Step 3): `Octave batch done: NH=449, NP=..., elapsed=...s.`
- The `--shambhala-timeout-s` appears in the job subprocess command (grep the log for `shambhala-timeout-s`).
- The first job reaches the S3 upload phase and exits with code 0.

6.4. If the first job succeeds, let the run continue to completion. Once all 54 jobs are done, verify that `shambhala_metrics.csv` is uploaded to S3 and that 1,188 output keys exist in `FL_batch_correction/exp/` matching the `shambhala_` pattern.

6.5. Download `shambhala_metrics.csv` and confirm the `r2_batch` values are in the expected range (well below the raw baseline of ~0.95).

---

### Completion Checklist

- [x] Step 1: `run_shambhala_parallel.py` — timeout propagation bug fixed, `--shambhala-timeout-s` CLI arg added
- [ ] Step 2: Dry-run verification that `--shambhala-timeout-s` appears in subprocess command *(run in pod)*
- [x] Step 3: Per-batch timing added to `parallel.py` `_harmonize_one_batch()`; `logging.basicConfig` added to `run_shambhala_job.py` so INFO lines are visible in job stdout
- [ ] Step 4: Single test job run on pod; NH and timing recorded in `k8s/README.md` timing table *(run in pod)*
- [x] Step 5: `k8s/README.md` updated with vCPU/workers guidance, two-level timeout explanation, and timing table
- [ ] Step 6: Full 54-job re-run completed; 1,188 S3 outputs verified; `shambhala_metrics.csv` downloaded and sanity-checked *(run in pod)*
