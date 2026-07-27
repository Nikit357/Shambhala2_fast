# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Package Is

The `shambhala/` directory is the installable Python library (`pip install -e .` from the repo root, `/home/jovyan/Projects/Shambhala2_fast`). It contains seven modules that together implement the Python side of the Shambhala2 harmonization pipeline. The Octave side lives in `../octave/`. The CLI entry point is `../run_shambhala.py`.

For full architecture, CLI usage, K8s operations, and benchmark integration, see `../CLAUDE.md`.

---

## Module Map

| Module | Sole responsibility |
|---|---|
| `octave_bridge.py` | Format TSV → pipe to Octave subprocess → parse stdout. No logic beyond I/O. |
| `parallel.py` | Split samples into batches → `ProcessPoolExecutor` → assemble → Q-rescale. |
| `q_rescale.py` | Compute Q reference statistics (`compute_q_statistics`); apply mean/std alignment (`rescale`). |
| `na_handling.py` | Drop or KNN-impute NaN genes before Octave; restore dropped genes as NaN after. |
| `io_utils.py` | Read/write `.csv`/`.tsv`/`.tsv.gz` from local paths or S3 URIs. |
| `progress_display.py` | `_ProgressRenderer` thread + `ProgressEvent` NamedTuple; ANSI TTY and plain-text modes. |
| `cublock_python.py` | `cublock_python` / `mod_pol_python` — NumPy+sklearn CuBlock for `--python-cublock`. Exploratory; deadlocks under `ProcessPoolExecutor` and diverges from Octave k-means. |

---

## Critical Invariants

### Data orientation
- **External API** (`io_utils`, `na_handling`, `q_rescale`): FL convention — rows = samples, columns = genes.
- **Octave bridge + `parallel.py` internals**: genes × samples (transposed). The transpose happens once at the top of `run_shambhala.py::main()`.
- Never transpose inside any of the six modules; each module's docstring states which convention it expects.

### Scale expectations (single most common bug source)
- `readExpressionData.m` applies `log2(x+1)` internally. Data piped to Octave via `octave_bridge.py` must be in **raw (non-log) scale**. If input is already log-transformed, exponentiate before calling `_format_pool_as_tsv`.
- Octave output is also raw scale. `q_rescale.rescale()` applies `log(x+1)` before the mean/std alignment, then exponentiates back.

### Why `ProcessPoolExecutor`, not `ThreadPoolExecutor`
Each worker runs an independent Octave subprocess. Threads share file descriptors; two threads piping to two Octave processes would corrupt each other's stdin. This constraint must be preserved.

### `_harmonize_one_batch` and `_harmonize_one_batch_with_progress` must stay module-level
`ProcessPoolExecutor` pickles the function before dispatch. Lambdas and closures are not picklable. Do not move these functions inside `harmonize_parallel`.

### `n_workers` cap
`harmonize_parallel` clamps `n_workers` to `min(n_workers, 10, n_samples)`. Never raise the hard cap of 10 without understanding file-descriptor limits on the target host.

### Q pseudocount
`compute_q_statistics` adds `q_pseudocount=1e-6` to all Q values before `log`. This prevents `-Inf` from zero-count RNA-seq genes. Passing `q_pseudocount=0` switches to warn-and-exclude mode. Do not silently drop the pseudocount.

---

## Running Tests

Tests live in `../tests/`. Always run from the repo root, not from inside this directory:

```bash
cd /home/jovyan/Projects/Shambhala2_fast
pytest tests/ -v                             # full suite (~109 tests)
pytest tests/test_octave_bridge.py -v        # single file
pytest tests/test_q_rescale.py::test_rescale_known_values -v  # single test
```

Octave-requiring tests auto-skip when Octave is unavailable (guarded by `@octave_required`). Pure-Python tests always run.

---

## Key Design Decisions Worth Preserving

- `octave_bridge.py` uses `run_octave_normalize_streaming` (not `run_octave_normalize`) in `parallel.py` so that per-sample `"Harmonizing sample X out of Y"` lines from Octave are intercepted and forwarded to `_ProgressRenderer` via a `multiprocessing.Manager().Queue()`.
- `disable_progress=True` in `harmonize_parallel` skips Manager creation entirely — use this in unit tests to avoid `multiprocessing.Manager` startup overhead.
- `_format_pool_as_tsv` writes `SYMBOL\t...` as the header. `readExpressionData.m` expects exactly that keyword in column 1.
- `rescale()` does an inner join on gene symbols between Octave output and `rm`/`rs`. Genes absent from Q statistics are dropped with a warning, not silently; downstream callers should not assume the output has the same gene set as the input.
