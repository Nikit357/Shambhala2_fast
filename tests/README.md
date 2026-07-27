# Tests — Shambhala_containerized

## Directory Layout

```
tests/
├── fixtures/                    # Input/output CSV fixtures (FL convention: rows=samples)
│   ├── input_10samples.csv      # 10 × 34 592
│   ├── P0_small.csv             # 39 × 34 592 (calibration reference)
│   ├── Q0_small.csv             # ~100 × 34 592 (definitive reference)
│   ├── expected_output.csv      # Baseline: n_workers=1, seed=42, no speed-up flags
│   ├── input_with_nas.csv       # 10 × 34 592 with ~5% NaN injected (seed=42)
│   └── generate_fixtures.py     # Regenerates all above from Shambhala2/ originals
├── speed_up_tests/              # Integration tests for optional speed-up flags
│   ├── conftest.py              # Shared: load_fixtures, run_pipeline, assert_mean_rel_diff
│   ├── test_approach_a.py       # Approach A: Python-side QN precomputation
│   ├── test_approach_b.py       # Approach B: skip P column polynomial fits (always active)
│   ├── test_approach_c.py       # Approach C: centroid P subsampling
│   ├── test_approach_d.py       # Approach D: Python-side QN (unit tests)
│   ├── test_approach_e.py       # Approach E: precomputed gene k-means clusters
│   ├── test_approach_g.py       # Approach G: full Python CuBlock (sklearn)
│   └── test_combinations.py    # Combined flag interactions
├── test_integration.py          # End-to-end baseline pipeline (no speed-up flags)
├── test_octave_bridge.py        # Octave subprocess I/O unit tests
├── test_q_rescale.py            # Q-rescaling unit tests
├── test_na_handling.py          # NaN handling unit tests
├── test_io_utils.py             # Local and S3 I/O unit tests
└── test_progress_display.py     # Progress renderer unit tests
```

## How to Run

Run all commands from the `Shambhala_containerized/` root directory.

```bash
# All tests — 1-hour timeout per test
pytest tests/ -v --timeout=3600

# Unit tests only — no Octave required, completes in < 30 s
pytest tests/ -v --ignore=tests/speed_up_tests/ --timeout=60

# Speed-up tests only — all, including timing tests (default)
pytest tests/speed_up_tests/ -v --timeout=3600

# Speed-up tests — skip timing tests for faster CI
pytest tests/speed_up_tests/ -v -m "not slow" --timeout=3600

# Speed-up tests — quick sanity check with 20-minute cap
pytest tests/speed_up_tests/ -v -m "not slow" --timeout=1200

# Single test file
pytest tests/test_q_rescale.py -v

# Single test
pytest tests/test_q_rescale.py::test_rescale_known_values -v
```

## Timeouts

Every speed-up test that calls Octave or Python CuBlock is decorated with
`@pytest.mark.timeout(3600)` (1 hour). A global fallback of 3600 s is set in
`pyproject.toml`. The `--timeout N` flag overrides both for the whole run.

```bash
pytest tests/speed_up_tests/ -v --timeout=1200   # 20-minute limit per test
pytest tests/speed_up_tests/ -v --timeout=0       # no limit (not recommended)
```

## Primary Integrity Metric

**`mean_rel_diff`** is the pass/fail criterion for all speed-up integrity tests.
It is the mean of `|result − expected| / |expected|` across all genes and samples,
with zero-protection on the denominator.

`max_abs_diff` is reported in every test output for diagnostic use only (to spot
outlier genes). It is not the basis for the assertion.

### Two-tier threshold

| Tier | `mean_rel_diff` | Test outcome | Interpretation |
|---|---|---|---|
| PASSED | < 1% | Green | Production-safe |
| WARNING | 1%–5% | Green + warning message | Exploratory only; not for benchmarking |
| FAILED | > 5% | Red | Unacceptable systematic distortion |

Approach A+synthetic uses wider tiers (warn > 10%, fail > 20%) because it is
an explicitly exploratory approximation with known large deviation.

### Threshold table by test

| Test | warn | fail | Expected `mean_rel_diff` |
|---|---|---|---|
| `test_a_qn_only` | 1% | 5% | ~0.34% → PASSED |
| `test_a_synthetic_cublock_p` | 10% | 20% | ~14.66% → WARNING (exploratory) |
| `test_b_output_identical_to_baseline` | — | `max_abs_diff < 1e-4` | ~0 |
| `test_e_output_within_tolerance` | 1% | 5% | TBD |
| `test_g_full_pipeline_within_tolerance` | 1% | 5% | TBD |
| `test_g_cublock_matches_octave` | 1% | 5% | TBD (log2 space) |

## Test Markers

| Marker | Applied to | Behavior |
|---|---|---|
| `@octave_required` | All tests needing Octave | Skipped silently when Octave is not found |
| `@pytest.mark.slow` | Timing tests (4 total) | Included by default; exclude with `-m "not slow"` |
| `@pytest.mark.timeout(N)` | All Octave/CuBlock tests | Hard wall-clock limit of N seconds |

## Fixture Regeneration

Requires `Shambhala2/` to be present as a sibling of `Shambhala_containerized/`:

```bash
python tests/generate_fixtures.py
```

## Approximate Wall Times (JupyterHub server, 10-sample fixture)

| Test group | Estimated time |
|---|---|
| Unit tests (no Octave) | < 30 s |
| `test_integration.py` | ~3–5 min |
| Speed-up tests, no timing tests | ~20–60 min |
| Speed-up tests, all (with timing) | ~60–120 min |
