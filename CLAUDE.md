# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

Pure-Python + Octave implementation of the **Shambhala2** gene expression harmonization algorithm (Borisov et al., 2021). Replaces the original R wrapper with a Python pipeline that calls Octave via stdin/stdout subprocess — no intermediate files, no R, no `rpy2`.

Two entry points:
- **`run_shambhala.py`** — single-run CLI (one input matrix → one harmonized output).
- **`harmonization_scripts/`** — cross-product benchmark: 3 imputations × 18 P/Q variants = 54 Shambhala runs → 1,512 S3 outputs (14 strategies × 2 post_rm each; J/K derived via Case B without Octave re-run).

This is method `20_shambhala` in the FL dissertation benchmark (`../harmonization-scripts/bench_shared.py`).

Repo root: `/home/jovyan/Projects/Shambhala2_fast`. The installed package name is
`shambhala-containerized` (see `pyproject.toml`), so older docs and log files call the repo
`Shambhala_containerized` — same thing.

`Calibration_datasets/` ships the four reference matrices used by the benchmark
(`P0_standard.csv.gz`, `Q0_standard.csv.gz`, `ANTE.csv.gz`, `OncoboxCancer.csv.gz`); the
S3 copies under `FL_batch_correction/calibration_datasets/` are the same files.

### Nested CLAUDE.md files — read the relevant one before editing that directory

| File | Covers |
|---|---|
| `shambhala/CLAUDE.md` | Module invariants: orientation, log scale, `ProcessPoolExecutor` picklability, Q pseudocount |
| `tests/CLAUDE.md` | Test-file map, Octave skip pattern, tolerance tiers |
| `tests/fixtures/CLAUDE.md` | Fixture shapes and provenance |
| `harmonization_scripts/CLAUDE.md` | 54-job cross-product benchmark, S3 layout, exit codes |
| `harmonization_scripts/k8s/README.md` | Pod lifecycle, rsync, troubleshooting |
| `reseach_implementation_plans/CLAUDE.md` | Index of every root-cause analysis and plan |

---

## Commands

```bash
# Install (required before running any scripts)
pip install -e .

# Configuration — required by everything under harmonization_scripts/.
# shambhala_bench_shared.py raises KeyError at import if this is unset (no default
# on purpose: a wrong-bucket default would write somewhere unintended). See .env.example.
export SHAMBHALA_S3_BUCKET=your-bucket

# Run all tests — Octave-requiring ones auto-skip if Octave unavailable
# Always run from the repo root (/home/jovyan/Projects/Shambhala2_fast); tests import the
# `shambhala` and `harmonization_scripts` packages relative to it.
pytest tests/ -v

# Unit tests only — no Octave, < 30 s
pytest tests/ -v --ignore=tests/speed_up_tests/ --timeout=60

# Run speed-up combination tests — use -s to see the summary table and expression comparisons
pytest tests/speed_up_tests/test_combinations.py -v -s --timeout=600

# Skip the timing-only tests
pytest tests/speed_up_tests/ -v -m "not slow" --timeout=3600

# Run a single test file
pytest tests/test_octave_bridge.py -v

# Run a single test
pytest tests/test_q_rescale.py::test_rescale_known_values -v

# Run CLI on toy 10-sample fixture
python run_shambhala.py \
    --input tests/fixtures/input_10samples.csv \
    --P tests/fixtures/P0_small.csv \
    --Q tests/fixtures/Q0_small.csv \
    --output /tmp/harmonized_test.csv \
    --octave-bin /opt/conda/bin/octave \
    --n-workers 1 \
    --random-seed 42 \
    --log-level DEBUG

# Cross-product benchmark — smoke test (1 variant, strict imputation, requires S3 + Octave)
python harmonization_scripts/run_shambhala_job.py \
    --imp strict \
    --method shambhala_P0std_Q0std \
    --n-shambhala-workers 2 \
    --octave-bin /usr/bin/octave \
    --random-seed 42 \
    --out-json /tmp/test_result.json

# Cross-product benchmark — full run (run from K8s pod, see harmonization_scripts/README.md)
nohup python harmonization_scripts/run_shambhala_parallel.py \
    --n-workers 3 --n-shambhala-workers 5 \
    --octave-bin /usr/bin/octave --skip-if-exists \
    --memory-limit-gb 10.0 --timeout-s 21600 --random-seed 42 \
    > /workspace/shambhala_run.log 2>&1 &

# Regenerate test fixtures (requires Shambhala2/ as sibling directory)
python tests/generate_fixtures.py
```

---

## Architecture

### Pipeline flow (`run_shambhala.py`)

```
read inputs (FL convention: samples×genes)
  → transpose to internal convention (genes×samples)
  → gene intersection (input ∩ P ∩ Q)
  → apply NA strategy on input.T (FL convention!) → clean_T = result.T
  → align p_T to clean_T.index (survivors after NA drop)
  → compute Q statistics from q_T.loc[clean_T.index].T (FL convention!)
  → harmonize_parallel(clean_T, p_T, ...)
       ProcessPoolExecutor → n_workers batches
       each batch: concat [batch | p_T] → Octave subprocess via stdin → parse stdout
       Q-rescaling applied once after all batches assembled
  → restore_na_genes(harmonized_T.T, na_mask) → result_df (samples×genes)
  → write output
```

### Data orientation convention — critical

| Context | Orientation | Shape |
|---|---|---|
| External API (`io_utils`, `na_handling`, `q_rescale`) | FL convention: rows=samples, columns=genes | (n_samples, n_genes) |
| Octave bridge (`octave_bridge.py`, `parallel.py`) | Internal: rows=genes, columns=samples | (n_genes, n_samples) |
| `readExpressionData.m` | Shambhala2 internal: rows=genes, first col=SYMBOL | same |

The transpose happens at the top of `main()`: `input_T = input_df.T`. **All `na_handling` and `q_rescale` calls must receive un-transposed data.**

### Module responsibilities

| Module | Role |
|---|---|
| `run_shambhala.py` | CLI entry point; orchestrates all steps |
| `shambhala/octave_bridge.py` | Only interface to Octave; formats TSV stdin, parses stdout; `run_octave_normalize_streaming` intercepts per-sample progress lines |
| `shambhala/parallel.py` | Splits samples into batches; `ProcessPoolExecutor` dispatch; Q-rescaling after assembly; wires up `_ProgressRenderer` |
| `shambhala/progress_display.py` | Live per-worker progress display (ANSI TTY mode and plain-text fallback); `ProgressEvent` NamedTuple; `_ProgressRenderer` thread |
| `shambhala/q_rescale.py` | Computes Q reference statistics; `rescale()` applies mean/std alignment |
| `shambhala/na_handling.py` | `apply_na_strategy` (drop/knn) and `restore_na_genes`; FL convention |
| `shambhala/io_utils.py` | Read/write local + S3, auto-detects .csv/.tsv/.tsv.gz |
| `shambhala/cublock_python.py` | Pure-NumPy/sklearn CuBlock (`cublock_python`, `mod_pol_python`) used by `--python-cublock`. Exploratory only — see xfail note below |
| `octave/Shambhala2_piped.m` | Modified Shambhala2.m: reads stdin, writes stdout, args injected via `--eval`; already emits `Harmonizing sample X out of Y` per sample |
| `octave/CuBlock.m` | Verbatim copy — zero changes |
| `octave/*.m` | Verbatim copies of `quantilenorm.m`, `kmeans.m`, `readExpressionData.m` |

**Octave script variants.** `octave_bridge.py` picks the `.m` file to `source()` from the
speed-up flags — there is no single entry script:

| Script | Selected when | Difference from `Shambhala2_piped.m` |
|---|---|---|
| `Shambhala2_piped.m` | default | — |
| `Shambhala2_piped_preqn.m` | `skip_qn=True` (`--precompute-qn-reference`) | `quantilenorm(EXP)` call removed — Python already did QN, so calling it again would double-normalize |
| `Shambhala2_piped_fixed.m` | `fixed_clusters` passed (`--precompute-cublock-clusters`) | calls `CuBlock_fixed(...)`; requires `FIXED_CLUSTERS` (G×1 int labels 1..k) in the `--eval` preamble |
| `CuBlock_fixed.m` | via `Shambhala2_piped_fixed.m` | skips k-means entirely, reusing precomputed gene clusters across all 30 repetitions |

`skip_qn` and `fixed_clusters` are checked in that order in both `run_octave_normalize` and
`run_octave_normalize_streaming`; keep the two selection blocks in sync when adding a variant.

### Progress display

`parallel.py` uses `_ProgressRenderer` from `progress_display.py` to show live per-worker progress during harmonization.

**TTY mode** (interactive terminal — `sys.stderr.isatty()` is True): renders N+4 fixed rows using ANSI escape codes, refreshed at up to 4 Hz:
```
[Shambhala] 181/500 samples | 5 workers
───────────────────────────────────────────────────────────────────────────────
 Worker 1    [██████████████░░░░░░░░░░░░░░░░]    36/100 ( 36%)  3.2 s/sample  ETA 3m 12s
 Worker 2    [████████████████░░░░░░░░░░░░░░]    40/100 ( 40%)  3.1 s/sample  ETA 2m 52s
 ...
───────────────────────────────────────────────────────────────────────────────
 Overall     [█████████████░░░░░░░░░░░░░░░░░]   181/500 ( 36%)  elapsed 9m 37s  ETA 17m 06s
```

**Non-TTY mode** (log file, K8s pod logs): emits plain timestamped lines to stderr instead:
```
[14:22:05][Worker 1]  20/100 (20%)  elapsed 1m 04s
[14:22:31][Worker 2]  20/100 (20%)  elapsed 1m 30s
```

The `harmonize_parallel()` function accepts `disable_progress=True` to skip renderer and Manager creation entirely (useful for unit tests).

### Octave bridge mechanics

Python builds an `--eval` preamble injecting `NH`, `NP`, `k`, and optional random seed, then calls `source('{octave_scripts_dir}/Shambhala2_piped.m')`. The full input matrix (batch samples + all P samples) is serialized as TSV and piped to Octave stdin. Octave writes `SYMBOL val1 val2 ...` lines to stdout.

**Why `ProcessPoolExecutor` not `ThreadPoolExecutor`:** Each worker runs an independent Octave subprocess. Threads would share file descriptors for stdin/stdout, corrupting the pipe.

### Cross-product benchmark (`harmonization_scripts/`)

**Key insight:** Shambhala normalizes each sample independently. Run once on the full S0 (no-removal) dataset; filter post-hoc to each strategy's sample set. This reduces 54 × 14 = 756 Octave runs to just 54.

| File | Role |
|---|---|
| `shambhala_bench_shared.py` | Constants (18 variants, 14 strategies, 3 imputations), S3 helpers, PCA R² metrics, `identify_outlier_batches`. Must NOT import `rpy2`. |
| `run_shambhala_job.py` | Worker: one `(imp × method)` → downloads S0, runs Shambhala, writes 28 S3 outputs + JSON sidecar. Exit codes: 0=ok, 1=harmonization failed, 2=low RAM, 3=S3 missing. |

`run_shambhala_job.py` routes each job through one of three cases at Step 0 (docstring at the top
of the file is authoritative):

- **Case A** — all 28 outputs already on S3 → early exit, nothing downloaded.
- **Case B** — `S0_no_removal__{imp}__{method}__post0` exists but some outputs are missing →
  `_run_reuse_path()` derives them by filtering that matrix. No Octave, no calibration download.
- **Case C** — S0 harmonized absent (or `--skip-if-exists` not passed) → full harmonization.

Cases A and B only apply with `--skip-if-exists`. `tests/test_s0_reuse.py` covers this routing.
| `run_shambhala_parallel.py` | Dispatcher: `ThreadPoolExecutor`, memory guard, failed-jobs log, consolidated `shambhala_metrics.csv`. |
| `derive_rare_batches_outputs.py` | Fast-path: derives `I_rare_batches_removed` outputs by filtering existing S0 harmonized matrices — no re-harmonization needed. Also uploads strategy I annotation files and re-consolidates `shambhala_metrics.csv`. |
| `k8s/shambhala-pod.yaml` | K8s pod manifest: `python:3.12-bookworm`, Octave only (no R), ~5-10 min startup, 16 vCPU / 128 GiB. |

S3 key pattern (identical to main benchmark): `FL_batch_correction/exp/{strat}__{imp}__{method}__post{0|1}.tsv.gz`

**NA strategy for benchmark:** Always `--na-strategy drop` (the default). S0 knn/softimpute datasets are already imputed — double-imputing distorts values.

---

## Key Constraints

- **CuBlock requires NaN-free input.** NA genes must be dropped or imputed before Octave.
- **`source()` vs `addpath()`:** `addpath` adds to the Octave function search path; `source()` resolves relative to cwd. Always use an absolute path in `source()` (`octave_bridge.py` line ~203).
- **`readExpressionData.m` applies `log2(x+1)` internally.** Data piped to Octave must be in raw (non-log) scale — exponentiate if input is already log-transformed. `q_rescale.rescale` expects raw-scale Octave output and applies `log(x+1)` before rescaling.
- **P samples are the quantile normalization anchor only.** They are concatenated with each batch for Octave and discarded from output. Align P to the same gene set as the cleaned input before concatenation.
- **`n_workers` capped at 10** (enforced in both CLI and `harmonize_parallel`).
- **`shambhala_bench_shared.py` must not import `rpy2`** — the parent `bench_shared.py` does, requiring a full R installation. All metric functions are copied verbatim.

---

## Tests

~109 test functions in three groups (see `tests/CLAUDE.md` for the per-file map):

| Group | Files | Octave needed |
|---|---|---|
| Unit | `test_io_utils`, `test_q_rescale`, `test_na_handling`, `test_octave_bridge`, `test_progress_display` | mostly no (mocked) |
| Integration | `test_integration.py` | yes |
| Benchmark logic | `test_s0_reuse.py` — `_run_reuse_path` and the Case A/B/C routing in `run_shambhala_job.py`, fully mocked S3 | no |
| Speed-up | `tests/speed_up_tests/` — one file per approach plus `test_combinations.py` | yes |

pytest config lives in `pyproject.toml`: global `timeout = 3600` (`timeout_method = "signal"`,
needs `pytest-timeout`) and a `slow` marker for the timing tests, which run by default.

### Fixtures

All fixtures in `tests/fixtures/` use FL convention (rows=samples). `tests/generate_fixtures.py` creates them from `../Shambhala2/Input.csv`, `P0.csv`, `Q0.csv`.

`expected_output.csv` is a **determinism regression baseline** produced by this pipeline itself
(`--random-seed 42 --n-workers 1`), not by the original Shambhala2.R.
`test_single_worker_matches_reference` checks against it with `atol=1e-4`. Regenerating fixtures
therefore rebaselines that test — do it only deliberately.

---

## Infrastructure (K8s pod)

Pod name: `shambhala-bench` in namespace `<your-namespace>`. Base image: `python:3.12-bookworm`. No R installed; Octave installed via apt (`octave octave-statistics`).

**K8s requirements vs local requirements.txt:** The pod's ConfigMap in `k8s/shambhala-pod.yaml` pins `pandas>=2.2,<3.0` and adds `psutil>=5.9` and `awscli>=1.44`. The local `requirements.txt` and `pyproject.toml` still say `pandas>=1.5` and lack `psutil`. These diverged intentionally (see `reseach_implementation_plans/pod_modules_error_reseach_plan_260515.md`). When editing pip deps for the pod, edit the ConfigMap in the YAML, not `requirements.txt`.

**Pod pip install command** (inside pod, after rsync):
```bash
cd /app/Shambhala_containerized
pip install -e .   # no --break-system-packages needed on python:3.12-bookworm
```

**`matlab` wrapper** installed in the pod at `/usr/local/bin/matlab`: wraps `octave --no-gui`, stripping MATLAB-only flags (`-nodesktop`, `-nosplash`, `-nodisplay`) that Shambhala2 passes internally via `system("matlab ...")`.

Full K8s operations guide: `harmonization_scripts/k8s/README.md`.

---

## Speed-up Approaches

Motivation: the NBGPL570 P reference (NP=250) is ~14–15× slower than P0std (NP=39) because both
QN and CuBlock scale with NP. Approach **B** (skipping wasted CuBlock polynomial fits for P
columns) is always on — it is baked into `CuBlock.m` and changes nothing numerically. Everything
else is opt-in.

### Flag → approach → plumbing

| CLI flag (`run_shambhala.py`) | Approach | Where it acts | `harmonize_parallel` param |
|---|---|---|---|
| `--precompute-qn-reference` | A | `main()`: builds the QN reference from P once via the `qnorm` library, normalizes `clean_T` in Python | `skip_qn=True` |
| `--synthetic-cublock-p` | D | `main()`: replaces all P columns with one synthetic column (mean of sorted P per gene rank) | (fewer P columns in `p_df`) |
| `--max-p-samples N` | C | `main()`: keeps the N P-samples closest to the P centroid | (smaller `p_df`) |
| `--precompute-cublock-clusters` | E | `precompute_clusters_from_p()` in `octave_bridge.py` (sklearn k-means on genes) | `fixed_clusters=<G×1 array>` |
| `--python-cublock` | G | replaces the Octave CuBlock call with `cublock_python` | `python_cublock=True` |

`--python-cublock` requires `--precompute-qn-reference` (guarded in `main()`). `--synthetic-cublock-p`
is only meaningful together with `--precompute-qn-reference`.

Results from the 2026-05-20 combination test run on the 10-sample fixture (`--timeout=600`, baseline ~48 s). Full details in `speed_up_test_results.md`.

| Combo | Flags | mean_rel_diff | Speedup | Status |
|---|---|---|---|---|
| baseline | — | 0.00% | 1.00x | PASS |
| A_only | `precompute_qn_reference=True` | 0.98% | ~1.03x | PASS |
| A_D | `precompute_qn_reference=True, synthetic_cublock_p=True` | 14.67% | ~5.3x | WARN (exploratory) |
| C_40 | `max_p_samples=40` | 0.00% | 1.00x | PASS (no effect at NP=39) |
| C_20 | `max_p_samples=20` | 2.45% | ~1.6x | WARN |
| E | `precompute_cublock_clusters=True` | 1.26% | ~10.6x | WARN |
| **A_E** | `precompute_qn_reference=True, precompute_cublock_clusters=True` | **1.39%** | **~8.6x** | **WARN — APPROVED** |
| G | `python_cublock=True` | — | — | XFAIL (deadlock) |
| A_D_E | `precompute_qn_reference=True, synthetic_cublock_p=True, precompute_cublock_clusters=True` | 9.00% | ~16.8x | PASS (exploratory threshold 20%) |
| A_D_G | `precompute_qn_reference=True, python_cublock=True, synthetic_cublock_p=True` | — | — | XFAIL (deadlock) |

**A_E is the approved speed-up for large long-running Shambhala benchmark runs** (e.g. NBGPL570 with NP=250). It combines Python QN reference caching and precomputed k-means clusters, achieving ~8.6x speedup with only 1.39% mean relative distortion — acceptable for downstream PCA R² and batch metrics.

G and A_D_G are permanently xfailed due to a `ProcessPoolExecutor` deadlock triggered by `python_cublock=True` (QueueFeederThread blocked on large-array result serialization). The Python CuBlock implementation also produces 226% mean_rel_diff because NumPy PCG64 and Octave's Mersenne Twister yield different k-means cluster sequences from the same seed.

**Speed-up test assertion rule:** Always use `assert_mean_rel_diff` from `conftest.py`. Never assert on `max_abs_diff` — Q-rescaled values can reach ~100 000 and max is dominated by a single outlier gene. Run with `-s` to see the per-combo summary table and expression comparison output.

---

## Incident History

All resolved issues are documented in `reseach_implementation_plans/CLAUDE.md`. Key entries:
- `failures_analysis_260514.md` — three bugs in `octave_bridge.py` / `run_shambhala.py` / `test_integration.py` (all fixed; all 30 tests pass).
- `pod_no_python_reseach_plan_260515.md` — Ubuntu 24.04 startup failures (resolved by switching to `python:3.12-bookworm`).
- `pod_modules_error_reseach_plan_260515.md` — pandas 1.5.x build failure on Python 3.12 (resolved by pinning pandas≥2.2 in K8s ConfigMap).
