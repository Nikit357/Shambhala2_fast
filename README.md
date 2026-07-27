# Shambhala2_fast

A pure-Python + Octave implementation of the **Shambhala2** gene expression harmonization algorithm. Replaces the original R wrapper with a self-contained Python pipeline that handles S3/local I/O, per-sample parallelization, and NA management, while keeping the Octave numerical core unchanged.

## Algorithm

Shambhala2 (Borisov et al., 2021, doi:10.1093/bioinformatics/btab105) harmonizes each input sample independently in three steps:

1. **Quantile normalization** — the input sample is jointly quantile-normalized together with all samples in the calibration reference dataset P (Bolstad et al., 2003, doi:10.1093/bioinformatics/19.2.185). P anchors the quantile distribution shape.
2. **CuBlock normalization** — genes are clustered into k groups by k-means; within each group a cubic polynomial is fitted to map the distribution to a reference shape. The process is repeated 30 times and averaged for robustness.
3. **Q-rescaling** — each gene in the output is linearly rescaled so that its mean and standard deviation match those of the same gene in the definitive reference dataset Q (computed in log space).

## Prerequisites

- Python 3.11+
- [Octave](https://octave.org/) — any version that supports `/dev/stdin` and `fprintf(1,...)`. Tested with Conda Octave at `/opt/conda/bin/octave`.
- Runtime Python packages: see `requirements.txt`

## Installation

```bash
pip install -e .
```

Or install runtime dependencies directly:

```bash
pip install -r requirements.txt
```

## Usage

```bash
python run_shambhala.py \
    --input  path/to/input.csv \
    --P      path/to/P_reference.csv \
    --output path/to/harmonized_output.csv \
    --Q      path/to/Q_reference.csv
```

### Real example

Runs out of the box on the toy fixture and the calibration references shipped in this
repository — no external data needed:

```bash
python run_shambhala.py \
    --input  tests/fixtures/input_10samples.csv \
    --P      Calibration_datasets/P0_standard.csv.gz \
    --Q      Calibration_datasets/Q0_standard.csv.gz \
    --output harmonized.csv \
    --n-workers 3 --random-seed 42 --log-level DEBUG
```

All file arguments accept local paths or S3 URIs (`s3://bucket/key`). Input CSV files must follow the **FL data convention**: rows = samples, columns = genes, first column = sample ID (used as index).

### Full CLI reference

```
usage: run_shambhala.py [-h]
                        --input PATH --P PATH --Q PATH --output PATH
                        [--k INT]
                        [--n-workers INT]
                        [--na-strategy {drop,knn}]
                        [--knn-k INT]
                        [--max-na-frac FLOAT]
                        [--q-pseudocount FLOAT]
                        [--octave-bin PATH]
                        [--timeout-s INT]
                        [--random-seed INT]
                        [--log-level {DEBUG,INFO,WARNING}]

required arguments:
  --input PATH          Input expression matrix. Rows=samples, columns=genes.
                        Formats: .csv, .tsv, .tsv.gz (auto-detected).
  --P PATH              Calibration reference P (quantile-normalization anchor).
  --Q PATH              Definitive reference Q (defines target expression shape).
  --output PATH         Destination for the harmonized output matrix (CSV).

optional arguments:
  --k INT               k-means gene clusters for CuBlock. Default: 5.
  --n-workers INT       Parallel worker processes (1–10). Default: 5.
  --na-strategy {drop,knn}
                        NaN handling: drop (default) removes genes with any NaN
                        (restored as NaN in output); knn imputes via KNNImputer.
  --knn-k INT           KNN neighbours for imputation. Default: 5.
  --max-na-frac FLOAT   Max NaN fraction before dropping a gene in knn mode. Default: 0.20.
  --q-pseudocount FLOAT Pseudocount added to Q before log (prevents log(0)). Default: 1e-6.
  --octave-bin PATH     Octave executable. Default: 'octave'. K8s: /opt/conda/bin/octave.
  --timeout-s INT       Per-batch Octave timeout in seconds. Default: 6000.
  --random-seed INT     Seed for reproducible k-means in CuBlock. Default: not set.
  --log-level           DEBUG / INFO / WARNING. Default: INFO.
```

### Examples

```bash
# Basic run with default settings
python run_shambhala.py \
    --input data/expression.csv \
    --P reference/P0.csv \
    --Q reference/Q0.csv \
    --output results/harmonized.csv

# Remote: explicit Octave path, 10 workers, fixed seed, S3 I/O
python run_shambhala.py \
    --input  s3://your-bucket/expression.tsv.gz \
    --P      s3://your-bucket/P0.csv \
    --Q      s3://your-bucket/Q0.csv \
    --output s3://your-bucket/harmonized.csv \
    --octave-bin /opt/conda/bin/octave \
    --n-workers 10 \
    --random-seed 42 \
    --log-level DEBUG
```

## Progress display

When `harmonize_parallel()` runs, a live progress display is written to `stderr`.

**Interactive terminal** (`stderr` is a TTY): renders fixed rows — one per worker — using ANSI escape codes, refreshed up to 4 Hz:

```
[Shambhala] 181/500 samples | 5 workers
───────────────────────────────────────────────────────────────────────────────
 Worker 1    [██████████████░░░░░░░░░░░░░░░░]    36/100 ( 36%)  3.2 s/sample  ETA 3m 12s
 Worker 2    [████████████████░░░░░░░░░░░░░░]    40/100 ( 40%)  3.1 s/sample  ETA 2m 52s
 Worker 3    [████████████████████░░░░░░░░░░]    50/100 ( 50%)  3.0 s/sample  ETA 2m 30s
───────────────────────────────────────────────────────────────────────────────
 Overall     [█████████████░░░░░░░░░░░░░░░░░]   181/500 ( 36%)  elapsed 9m 37s  ETA 17m 06s
```

**Log file / K8s pod** (`stderr` is not a TTY): emits plain timestamped lines at most ~20 times per worker:

```
[14:22:05][Worker 1]  20/100 (20%)  elapsed 1m 04s
[14:22:31][Worker 2]  20/100 (20%)  elapsed 1m 30s
```

Per-sample granularity is possible because `Shambhala2_piped.m` already emits `Harmonizing sample X out of Y` for each sample; the Python bridge intercepts these lines from the Octave subprocess stdout without any Octave changes.

To suppress the progress display (e.g., in scripts that already capture stderr), call:
```python
harmonize_parallel(..., disable_progress=True)
```

## Calibration datasets

`Calibration_datasets/` ships four already-published reference matrices. Three are P
references (the quantile-normalization anchor); one is the Q reference (the definitive target
distribution). All use FL convention: rows = samples, columns = genes.

| File | Role | Samples | Genes | Description |
|---|---|---|---|---|
| `P0_standard.csv.gz` | P | 39 | 11,768 | Standard Shambhala2 calibration set (GEO `GSM…`) |
| `ANTE.csv.gz` | P | 202 | 36,596 | Normal-tissue reference, tissue-labelled |
| `OncoboxCancer.csv.gz` | P | 779 | 36,596 | Cancer-sample reference |
| `Q0_standard.csv.gz` | Q | 100 | 11,887 | Standard definitive reference (GTEx) |

Gene sets differ between files; the pipeline intersects input ∩ P ∩ Q automatically, so any
P/Q pair works.

### Standard calibration (P0 + Q0)

```bash
python run_shambhala.py \
    --input  tests/fixtures/input_10samples.csv \
    --P      Calibration_datasets/P0_standard.csv.gz \
    --Q      Calibration_datasets/Q0_standard.csv.gz \
    --output harmonized_P0std_Q0std.csv \
    --n-workers 3 --random-seed 42
```

### Normal-tissue anchor (ANTE + Q0)

A larger, more diverse P reference. Slower — quantile normalization and CuBlock both scale
with the number of P samples.

```bash
python run_shambhala.py \
    --input  tests/fixtures/input_10samples.csv \
    --P      Calibration_datasets/ANTE.csv.gz \
    --Q      Calibration_datasets/Q0_standard.csv.gz \
    --output harmonized_ANTE_Q0std.csv \
    --n-workers 3 --random-seed 42
```

### Cancer anchor (OncoboxCancer + Q0)

The largest P reference at 779 samples. Use `--precompute-qn-reference` to keep runtime
manageable.

```bash
python run_shambhala.py \
    --input  tests/fixtures/input_10samples.csv \
    --P      Calibration_datasets/OncoboxCancer.csv.gz \
    --Q      Calibration_datasets/Q0_standard.csv.gz \
    --output harmonized_Oncobox_Q0std.csv \
    --n-workers 3 --random-seed 42 \
    --precompute-qn-reference
```

### Loading a reference in Python

```python
from shambhala.io_utils import read_expression

p_reference = read_expression("Calibration_datasets/P0_standard.csv.gz")
print(p_reference.shape)   # (39, 11768) — rows = samples, columns = genes
```

### Original Shambhala2 references

The datasets from the Shambhala2 paper are on Zenodo: https://zenodo.org/record/6415067.
They use the Shambhala2-internal convention (rows = genes). The fixture files in
`tests/fixtures/` are transposed versions in FL convention (rows = samples) — the orientation
this pipeline expects.

## Configuration

The single-run CLI (`run_shambhala.py`) needs no configuration. The benchmark harness under
`harmonization_scripts/` reads its S3 bucket from the environment and raises at import time if
it is unset:

```bash
export SHAMBHALA_S3_BUCKET=your-bucket
```

See `.env.example`. There is deliberately no default — a wrong-bucket default would fail
silently or write somewhere unintended.

## Testing

Generate fixtures once (requires Shambhala2/ to be present as a sibling directory):

```bash
python tests/generate_fixtures.py
```

Run all tests:

```bash
pytest tests/ -v
```

Tests that require Octave are automatically skipped when Octave is not installed. All other tests (unit tests for IO, NA handling, Q-rescaling, and mocked Octave bridge) always run.

### Test matrix

#### `test_integration.py` — end-to-end tests (require Octave)

| Test | Description | Pass criterion |
|---|---|---|
| `test_single_worker_matches_reference` | n_workers=1 on 10-sample fixture | Output matches `expected_output.csv` within atol=1e-4 |
| `test_parallel_matches_single` | n_workers=3 and n_workers=5 | Per-gene row means match n_workers=1 within atol=0.05 |
| `test_local_io_roundtrip` | Write to temp file; read back | Shape and values identical |
| `test_gene_intersection_is_applied` | P has fewer genes than input | Output gene count = intersection size |
| `test_na_drop_full_pipeline` | Input with 5% NaN, strategy=drop | NaN genes remain NaN in output; others match expected |
| `test_na_knn_full_pipeline` | Input with 5% NaN, strategy=knn | Output has no NaN; shape matches input |

#### `test_na_handling.py` — unit tests

| Test | Description |
|---|---|
| `test_detect_na_genes_correct` | 3 known NaN genes → returns exactly those 3 |
| `test_drop_removes_na_genes` | strategy=drop → returned df has zero NaN |
| `test_knn_fills_nas` | strategy=knn → returned df has zero NaN |
| `test_restore_preserves_original_order` | Column order matches original after restore |
| `test_restore_inserts_nan_correct_columns` | Dropped gene columns are NaN in restored output |
| `test_all_na_gene_handled_gracefully` | 100% NaN gene with knn strategy → dropped, no exception |

#### `test_octave_bridge.py` — unit tests

| Test | Requires Octave | Description |
|---|---|---|
| `test_bridge_returns_correct_shape` | Yes | 5-sample batch + P → output (n_genes, 5) |
| `test_bridge_symbol_index_correct` | Yes | Row index matches gene symbols from input |
| `test_bridge_raises_on_nonzero_exit` | No | Non-existent binary → OctaveExecutionError |
| `test_bridge_error_message_contains_stderr` | No (mocked) | Exception message includes stderr |
| `test_bridge_raises_on_unparseable_stdout` | No (mocked) | Garbage stdout → OctaveOutputError |
| `test_tsv_formatting_round_trip` | No | _format_pool_as_tsv → read back → values match |
| `test_random_seed_injected_in_preamble` | No (mocked) | --eval string contains rand('state', 42) |

#### `test_q_rescale.py` — unit tests

| Test | Description |
|---|---|
| `test_rescale_known_values` | Toy input with hand-computed output → matches within 1e-10 |
| `test_compute_q_statistics_shape` | Returns two Series indexed by gene symbol |
| `test_pseudocount_prevents_log_zero` | Zero Q value + default pseudocount → no warning |
| `test_pseudocount_zero_excludes_with_warning` | Zero Q value + pseudocount=0 → warning + exclusion |
| `test_rescale_raises_on_negative_q` | Negative Q value → ValueError |
| `test_rescale_gene_intersection` | 5 genes in input, 4 in Q → 4 genes in output |

#### `test_io_utils.py` — unit tests

| Test | Description |
|---|---|
| `test_read_local_csv` | Read fixture CSV → correct shape and dtype |
| `test_read_local_tsv_gz` | Read gzipped TSV → same result as uncompressed |
| `test_write_and_read_roundtrip` | Write df to temp file; read back → identical |
| `test_parse_s3_uri` | 's3://b/k/f.csv' → ('b', 'k/f.csv') |
| `test_s3_read_mocked` | Mock boto3 get_object; assert correct bucket/key |

---

## Speed-Up Options

The NBGPL570 P reference (NP=250) is ~14–15× slower than P0std (NP=39) because both QN and CuBlock scale with NP. Several optional flags reduce runtime with controlled accuracy trade-offs. **Approach B is always active** (baked into `CuBlock.m`).

### Flag overview

| Flag | Approach | Default | Integrity ⭐ | Description |
|---|---|---|---|---|
| *(none)* | **B** | Always on | ⭐⭐⭐⭐⭐ | Skip wasted CuBlock polynomial fits for P columns. Zero numerical change. |
| `--precompute-qn-reference` | **A+D** | Off | ⭐⭐⭐⭐ | Compute QN reference from P once in Python; apply per-sample via `qnorm`. Minor numerical deviation (< 0.05). |
| `--synthetic-cublock-p` | **A** extra | Off | ⭐⭐⭐ | Replace P columns in Octave with a single synthetic column (mean of QN reference). Only useful with `--precompute-qn-reference`. |
| `--max-p-samples N` | **C** | Off | ⭐⭐⭐⭐ | Subsample P to N samples closest to P centroid. Reduces QN anchor diversity. |
| `--precompute-cublock-clusters` | **E** | Off | ⭐⭐⭐ | Compute gene k-means clusters from P once using sklearn; inject into `CuBlock_fixed.m`. Clusters differ from per-sample Octave k-means. |
| `--python-cublock` | **G** | Off | ⭐⭐⭐ | Replace Octave CuBlock entirely with Python (sklearn + numpy). Fastest option; k-means init differs from Octave. |

⭐⭐⭐⭐⭐ = Identical to baseline · ⭐⭐⭐⭐ = Negligible numerical change · ⭐⭐⭐ = Acceptable deviation (< 0.10 max abs diff)

### Recommended combinations

| Goal | Flags | Expected speedup (NBGPL570) |
|---|---|---|
| Maximum fidelity | *(none)* | 1× (Approach B already active) |
| Fast + faithful | `--precompute-qn-reference` | 2–4× |
| Fast + faithful + small P | `--precompute-qn-reference --max-p-samples 40` | 4–8× |
| Maximum speed | `--precompute-qn-reference --synthetic-cublock-p --python-cublock` | 10–20× |
| Full Python (no Octave) | `--python-cublock` | Variable; validates the Python stack |

### Validation

Run the combination test suite on the K8s pod after enabling new flags:

```bash
# All speed-up tests — timing tests included (default)
pytest tests/speed_up_tests/ -v --timeout=3600

# Skip timing tests for faster CI
pytest tests/speed_up_tests/ -v -m "not slow" --timeout=3600

# Quick sanity check with 20-minute cap
pytest tests/speed_up_tests/ -v -m "not slow" --timeout=1200
```

**Primary integrity metric: `mean_rel_diff`** — mean of `|result − expected| / |expected|` across all genes and samples.

| Tier | `mean_rel_diff` | Test outcome | Interpretation |
|---|---|---|---|
| PASSED | < 1% | Green | Production-safe |
| WARNING | 1%–5% | Green + warning | Exploratory only; not for benchmarking |
| FAILED | > 5% | Red | Unacceptable systematic distortion |

Approach A+synthetic uses wider tiers (warn > 10%, fail > 20%) because it is an explicitly exploratory approximation.

`max_abs_diff` is reported alongside `mean_rel_diff` in every test for diagnostic purposes (spotting outlier genes) but is **not** the basis for the pass/fail assertion. In raw expression scale (post Q-rescaling), values reach tens of thousands, making absolute thresholds misleading.

---

## License

MIT — see `LICENSE`.

The Octave sources under `octave/` are vendored third-party or derived work and keep their own
upstream terms; the calibration datasets are previously published data redistributed here for
reproducibility. See `THIRD_PARTY_LICENSES.md` for the details and attributions.

## Citation

This is a reimplementation of Shambhala2:

> Borisov N, et al. Shambhala-2: a protocol for uniformly shaped harmonization of gene
> expression profiles of various formats. *Bioinformatics* (2021).
> doi:[10.1093/bioinformatics/btab105](https://doi.org/10.1093/bioinformatics/btab105)

The CuBlock normalization step is from:

> Junet V, Farrés J, Mas JM, Daura X. CuBlock: a cross-platform normalization method for
> gene-expression microarrays. *Bioinformatics* (2021).
