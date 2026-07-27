# CLAUDE.md — Research Implementation Plans Index

This directory holds root-cause analyses and implementation plans for the `Shambhala_containerized/` pipeline.
All plans are referenced relative to `Shambhala_containerized/`.

---

## 1. `failures_analysis_260514.md` — Test Suite Bug Fixes

**Date:** 2026-05-14  
**Status:** Complete — all 30 tests pass.

### What happened
`pytest tests/ -v` produced 8 failed / 22 passed. Three independent bugs.

### Bugs and fixes

| Bug | Location | Root cause | Fix |
|---|---|---|---|
| Bug 1 | `shambhala/octave_bridge.py:203` | `source('Shambhala2_piped.m')` uses a relative path; Octave resolves it from cwd (project root), not from the `octave/` dir. | Change to `source('{octave_scripts_dir}/Shambhala2_piped.m')` — absolute path. |
| Bug 2 | `run_shambhala.py` Steps 5 & 8; `tests/test_integration.py::_run_pipeline` | `apply_na_strategy` and `restore_na_genes` received genes×samples (transposed) instead of FL convention (samples×genes). Drop strategy then treated sample names as genes, dropping all 10 samples → 0-column input → `ValueError` in `harmonize_parallel`. | Pass `input_T.T` to `apply_na_strategy`; re-transpose result to `clean_T`. Pass `harmonized_T.T` to `restore_na_genes`. |
| Bug 3 | `run_shambhala.py` Step 6; `tests/test_integration.py::_run_pipeline` | `compute_q_statistics` received genes×samples. `gene_symbols = q_T.columns` returned 100 GTEX sample IDs instead of 8174 gene names. Q statistics indexed by sample names → empty intersection with Octave output → all-empty harmonized matrices. | Pass `q_T.loc[clean_T.index].T` (FL convention, filtered to NA-surviving genes). |

All three bugs were also present in `tests/test_integration.py::_run_pipeline` and fixed in the same pass.

### Result
30/30 tests pass. `test_single_worker_matches_reference` confirms output matches original Shambhala2.R within `atol=1e-4`.

---

## 2. `pod_no_python_reseach_plan_260515.md` — Pod Startup: pip/python Missing (Attempt 1)

**Date:** 2026-05-15  
**Status:** Applied — superceded by attempt 2 (see §3) and python image switch (see §5).

### What happened
First pod deployment: `python` and `pip3` not found despite `[startup] Python packages installed` appearing in logs. Startup looked successful but pod was non-functional.

### Root causes

1. **`liboctave-dev` removed from Ubuntu 24.04.** The apt-get command included `liboctave-dev`, which does not exist in Ubuntu 24.04. When apt encounters an unresolvable package it aborts the entire command without installing anything — including `python3-pip` and `python-is-python3`.

2. **No `set -e` in startup script.** Each `echo "[startup] ... installed"` line ran unconditionally, reporting false success after failed install commands.

### Fix applied
- Removed `liboctave-dev` from the step-2 apt-get command.
- Added `set -e` as the first line of the startup bash script.
- Added `octave --version` and `python3 --version && pip3 --version` as verification lines after install.

**Note:** `liboctave-dev` was never needed — Shambhala calls Octave only via stdin/stdout subprocess; it never links against Octave C++ headers.

---

## 3. `pod_modules_error_reseach_plan_260515.md` — Pod Startup: pandas Build Failure (Attempt 2)

**Date:** 2026-05-15  
**Status:** Applied — requirements.txt updated; superceded by python image switch (see §5).

### What happened
After fixing attempt 1, pip install failed building pandas 1.5.3 from source with `ModuleNotFoundError: No module named 'pkg_resources'`.

### Root causes

1. **pandas 1.5.3 has no Python 3.12 wheel.** pandas 1.5.x predates Python 3.12. pip falls back to downloading the source tarball and building it. The build invokes `setup.py` which imports `pkg_resources` from `setuptools`, but pip's build-isolation sandbox does not expose it correctly on Ubuntu 24.04's pip 24.

2. **`awscli>=1.29` (pip) conflicts with boto3/botocore.** awscli 1.x pins boto3/botocore to specific older versions, conflicting with the unconstrained boto3/botocore bounds in requirements.txt.

### Fix applied
- `pandas>=1.5.3,<2.0` → `pandas>=2.2,<3.0` — pandas 2.2 has a pre-built cp312 wheel; no source build.
- `awscli` removed from pip requirements, moved to apt-get in the startup script (Ubuntu 24.04 apt version is pre-pinned to compatible boto3/botocore).

### Compatibility check
All pandas APIs used by the pipeline (`read_csv`, `concat`, `loc`, `.T`, `.index.intersection`, `.fillna`) are stable across 1.5→2.x. No code changes required.

---

## 4. `cross_product_implementation_plan_260515.md` — Cross-Product Benchmark (54 runs → 1,188 S3 outputs)

**Date:** 2026-05-15  
**Status:** Code complete (Phases 0–5). Phase 6 (pod validation, S3 run) pending cluster access.

### Goal
Extend `run_shambhala.py` into a full benchmark matching the `harmonization-scripts/` framework. Key insight: Shambhala normalizes each sample independently, so 11 batch-removal strategies do not require 11 separate Shambhala runs — filter post-hoc.

### Dimensions

| Axis | Count |
|---|---|
| Imputation (strict / knn / softimpute) | 3 |
| Shambhala variants (9 P × 2 Q combinations) | 18 |
| Batch-removal strategies | 11 |
| Post-removal (post0 / post1) | 2 |
| **Total S3 outputs** | **1,188** |
| **Total Shambhala runs** | **54** |

### Files created

| File | Role |
|---|---|
| `harmonization_scripts/shambhala_bench_shared.py` | Constants, 18-variant registry, S3 key functions, post-removal, PCA R² metric. Must not import rpy2. |
| `harmonization_scripts/run_shambhala_job.py` | Worker: one (imp × method) → 22 S3 outputs + JSON sidecar. Exit codes: 0=ok, 1=harmonization failed, 2=low RAM, 3=S3 missing. |
| `harmonization_scripts/run_shambhala_parallel.py` | Dispatcher: `ThreadPoolExecutor`, memory guard, failed-jobs log, consolidated `shambhala_metrics.csv`. |
| `harmonization_scripts/k8s/shambhala-pod.yaml` | K8s pod manifest. |
| `harmonization_scripts/k8s/README.md` | Pod launch, SSH, rsync, teardown guide. |
| `harmonization_scripts/README.md` | Full CLI reference, timing estimates, NA strategy note. |

### S3 output key pattern
Identical to main benchmark: `FL_batch_correction/exp/{strat}__{imp}__{method}__post{0|1}.tsv.gz`
Fully compatible with `harmonization-metrics/` pipeline without modification.

### Phase completion

| Phase | Description | Status |
|---|---|---|
| 0 | Fix 3 bugs in `run_shambhala.py` / `octave_bridge.py` | ✅ Done |
| 1 | `shambhala_bench_shared.py` | ✅ Done |
| 2 | `run_shambhala_job.py` | ✅ Done (S3 smoke test pending pod) |
| 3 | `run_shambhala_parallel.py` | ✅ Done (S3 test pending pod) |
| 4 | K8s pod manifest + README | ✅ Done (pod deploy pending) |
| 5 | `harmonization_scripts/README.md` | ✅ Done |
| 6 | End-to-end validation on pod (1,188 S3 outputs) | ⏳ Pending cluster access |

---

## 5. `../python_docker_plan.md` — Switch Base Image to python:3.12-bookworm

**Date:** 2026-05-16  
**Status:** Applied to `harmonization_scripts/k8s/shambhala-pod.yaml`.

### What happened
The pod crashed on Ubuntu 24.04 across two successive attempts (§2 and §3 above). The root of both failures was that Ubuntu 24.04's packaging environment conflicted with the Python toolchain setup:
- `liboctave-dev` absent from Ubuntu 24.04 broke the entire apt step.
- Ubuntu's PEP 668 protection required `--break-system-packages`.
- pandas 1.5.x had no cp312 wheel, forcing a broken source build.

Switching to the official Python image eliminates all of these problems.

### Selected image: `python:3.12-bookworm`

Built on `buildpack-deps:bookworm` (Debian 12). Pre-installs: Python 3.12, pip, build-essential, curl, wget, gnupg, ca-certificates. Octave packages available via apt (same Debian 12 repos). Not alpine (Octave unavailable/broken musl). Not slim (would need to reinstall build-essential etc.).

### Changes applied

| Location | Before | After |
|---|---|---|
| `image:` | `ubuntu:24.04` | `python:3.12-bookworm` |
| Section-1 apt | `openssh-server rsync curl wget gnupg ca-certificates tmux htop` | `openssh-server rsync tmux htop` |
| Section-2 apt | `build-essential python3 python3-pip python3-dev python-is-python3 octave octave-statistics awscli` | `octave octave-statistics` |
| pip install | `pip3 install --break-system-packages --no-cache-dir` | `pip install --no-cache-dir` |
| verification | `python3 --version && pip3 --version` | `python --version && pip --version` |
| ConfigMap requirements.txt | no awscli entry | `awscli>=2.13` added (pip install; 2.x release) |
| `CLAUDE.md` pod table | `ubuntu:24.04` | `python:3.12-bookworm` |

### Result
- `apt-get install` reduced from ~12 to 6 packages.
- `--break-system-packages` flag removed (official Python image pip is not PEP 668 managed).
- pandas 2.2 (already in requirements from §3 fix) installs via pre-built cp312 wheel without issue.
- No functional changes to Octave setup, SSH, or rsync.

---

## 6. `shambhala_failures_research_implementation_260516.md` — Pod Run Failures (2026-05-16)

**Date:** 2026-05-16  
**Status:** Root causes identified; fixes not yet applied.

### Issues found

| Issue | Severity | Root Cause |
|---|---|---|
| GTExAffy P dataset: 0 gene intersection | Critical | CSV stored in genes×samples (Shambhala2 convention), not FL convention (samples×genes) |
| No guard against empty intersection | Critical | Missing Python check; Octave crashes with `r_sort(5): out of bound 0` |
| P/Q datasets not checked for NaN | Moderate | No assertion after `download_calib_df`; hidden NaN would silently corrupt Octave |
| 31–34% gene loss for Zenodo P/Q refs | Moderate | S0 strict (3447 genes) ≠ Zenodo gene universe; QNBKass avoids this for ANTE |
| DtypeWarning in annotation loading | Minor | `low_memory=False` missing in `pd.read_csv` |

### Key finding on knn strategy

`--na-strategy knn` does NOT help with the observed crashes. Strict S0 has 0 NaN in expression (confirmed: "0 genes dropped"). The crashes are caused by wrong GTExAffy orientation and missing Python guard. However, checking P and Q for NaN separately (they bypass the existing NA strategy) is the correct defensive fix.

### Fixes required (in priority order)

1. Re-prepare `GTExAffymetrix.csv` in FL convention and re-upload to S3.
2. Add empty-intersection guard in `run_shambhala_job.py` and `run_shambhala.py`.
3. Add defensive `n_probes == 0` check in `kmeans.m`.
4. Add NaN assertions for P and Q after `download_calib_df`.
5. Add `low_memory=False` to `download_ann_from_s3`.
