# Implementation Plan: `Shambhala_containerized`

**Date:** 2026-05-14  
**Author:** Daniil Nikitin / Claude Code  
**Status:** Plan — approved by Daniil; ready for implementation  
**Context:** See `Shambhala_research_260513.md` for background on the Shambhala2 algorithm, reference datasets (P/Q), and the original containerization design study.

---

## Overview

`Shambhala_containerized` is a fully self-contained, production-grade harmonization pipeline built around the Shambhala2 algorithm. It replaces the original R-centric wrapper with a pure-Python entry point that handles S3/local I/O, per-sample parallelization, and NA management, while keeping the Octave numerical core as close to the original as physically possible (exactly 3 line-level edits to `Shambhala2.m`; all other `.m` files copied verbatim).

The pipeline is designed to slot directly into the FL dissertation benchmark via its standardized samples × genes I/O convention.

---

## 0. Guiding Principles

| Requirement | Constraint | Architectural Decision |
|---|---|---|
| No intermediate files | Data must flow in-memory across the full pipeline | Python → Octave via `subprocess` stdin/stdout |
| Minimal Octave changes | Preserve algorithmic integrity, prevent AI-introduced deviations | Exactly 3 targeted line edits to `Shambhala2.m`; zero changes to all other `.m` files |
| Parallelization | Each Shambhala sample is independent; `rpy2` is not thread-safe | `ProcessPoolExecutor`; each worker is an isolated subprocess calling Octave |
| Eliminate R from the pipeline | The R bridge only does simple matrix ops (merge, log, rescale) | Replace with pure Python/pandas/numpy; removes `rpy2` dependency entirely, enabling clean multi-process execution |
| FL data convention | FL project uses rows = samples, columns = genes | Transpose to genes × samples on entry to Octave, transpose back on exit |
| S3 + local I/O | Input/output can be either local paths or `s3://` URIs | Unified reader/writer detecting `s3://` prefix |
| NA handling | CuBlock docstring: "The micro array cannot contain NaN values" | Two configurable strategies: `drop` (default) or `knn` (via `--na-strategy knn`); applied before Octave, restored/preserved after |

**Note on R elimination:** The original `Shambhala2/Shambhala2.R` and its `Shambhala2_flex` entry point used by the benchmark pipeline via `rpy2` in `bench_shared.py` are **not touched**. The new `Shambhala_containerized/` is a fully independent directory. The only reason R existed in the original pipeline was to merge DataFrames and apply the final Q-rescaling — both are trivial pandas/numpy operations.

---

## 1. Directory Structure

```
shambhala_adoption/
├── Shambhala2/                        # Original, untouched
└── Shambhala_containerized/           # New — this project
    │
    ├── run_shambhala.py               # Single CLI entry point
    ├── README.md                      # Detailed repository description and full usage guide
    ├── pyproject.toml                 # Package metadata and build config
    ├── requirements.txt               # Pinned runtime dependencies
    ├── setup.py                       # Compatibility shim for editable installs
    │
    ├── shambhala/                     # Core library package
    │   ├── __init__.py
    │   ├── octave_bridge.py           # Python → Octave subprocess interface
    │   ├── q_rescale.py               # Q-rescaling step (Python replacement of R logic)
    │   ├── na_handling.py             # NA detection, drop strategy, KNN imputation
    │   ├── parallel.py                # Parallel sample-batch dispatcher
    │   └── io_utils.py                # Unified S3 / local file I/O
    │
    ├── octave/                        # Octave scripts
    │   ├── Shambhala2_piped.m         # Modified Shambhala2.m — exactly 3 line-level changes
    │   ├── CuBlock.m                  # Verbatim copy from Shambhala2/ — ZERO changes
    │   ├── quantilenorm.m             # Verbatim copy — ZERO changes
    │   ├── kmeans.m                   # Verbatim copy — ZERO changes
    │   └── readExpressionData.m       # Verbatim copy — ZERO changes
    │
    └── tests/
        ├── generate_fixtures.py       # Script to produce all fixture files from Shambhala2/ data
        ├── test_integration.py        # Full pipeline end-to-end tests
        ├── test_na_handling.py        # Unit tests for NA module
        ├── test_octave_bridge.py      # Unit tests for the Octave bridge
        ├── test_q_rescale.py          # Unit tests for Q-rescaling math
        ├── test_io_utils.py           # Unit tests for S3 + local I/O
        └── fixtures/
            ├── input_10samples.csv    # Toy input: 10 samples × ~2000 genes (FL convention)
            ├── expected_output.csv    # Ground-truth output for the toy input (FL convention)
            ├── input_with_nas.csv     # Same toy input with ~5% synthetic NaN values
            ├── P0_small.csv           # P calibration reference (FL convention); transposed from Shambhala2/P0.csv
            └── Q0_small.csv           # Q definitive reference (FL convention); transposed from Shambhala2/Q0.csv
```

All fixture CSV files use the **FL data convention** (rows = samples, columns = genes). They are derived from the existing `Shambhala2/Input.csv`, `Output.csv`, `P0.csv`, and `Q0.csv` by transposing.

**README.md must include:**
- Algorithm description (three-step Shambhala2 process)
- Installation instructions (`pip install -e .`)
- Full CLI usage guide with all flags and examples
- Description of all test cases (see §5) so that Daniil can verify coverage before running

---

## 2. Changes to Octave Scripts

### 2.1 Summary

| File | Action | Changes |
|---|---|---|
| `Shambhala2.m` | Modified → `Shambhala2_piped.m` | Exactly 3 targeted edits (described below) |
| `CuBlock.m` | Verbatim copy | None |
| `quantilenorm.m` | Verbatim copy | None |
| `kmeans.m` | Verbatim copy | None |
| `readExpressionData.m` | Verbatim copy | None — works with `/dev/stdin` already |

### 2.2 Three Edits to `Shambhala2_piped.m`

#### Change 1 — Remove `args.txt` reading block

**Variable definitions:**
- `NH` — number of input samples in this batch (the batch that this Octave subprocess will normalize; for parallelization, one subprocess handles a subset of the full dataset).
- `NP` — number of P calibration reference samples.
- `k` — number of k-means clusters for the CuBlock algorithm (default 5).

The Python subprocess passes `NH`, `NP`, and `k` directly by prepending them in the `--eval` string before sourcing the script. The entire 5-line block that reads `args.txt` is deleted.

```matlab
% DELETE THESE 5 LINES (NH, NP, k arrive from --eval preamble):
fileID = fopen('args.txt','r');
sizeA = [3 1];
A = fscanf(fileID,'%f',sizeA);
fclose(fileID);
NH = A(1); NP = A(2); k = A(3);
```

No replacement text is needed. `NH`, `NP`, and `k` will already exist in the Octave workspace when the script is sourced.

#### Change 2 — Read pool data from stdin instead of `P_prim.txt`

```matlab
% ORIGINAL:
inData=readExpressionData("P_prim.txt",'log2');

% REPLACE WITH:
inData=readExpressionData('/dev/stdin','log2');
```

`readExpressionData.m` already calls `fopen(filename,'r')`. Passing `/dev/stdin` causes it to open standard input as a file descriptor, which `textscan` reads sequentially on Linux. **Zero changes are required to `readExpressionData.m` itself.**

#### Change 3 — Write output to stdout instead of `Cu_bis.txt`

```matlab
% ORIGINAL (the entire file-write block, ~7 lines):
outFileName = 'Cu_bis.txt';
outFile=fopen(outFileName,'w');
nG=size(OUT,1);
nS=size(OUT,2);
for i=1:nG
    fprintf(outFile,'%s',SYMBOL{i,1});
    for j=1:nS
        fprintf(outFile,' %f',OUT(i,j));
    end
    fprintf(outFile, '\n');
end
fclose(outFile);

% REPLACE WITH (identical loop body, file descriptor 1 = stdout):
nG=size(OUT,1);
nS=size(OUT,2);
for i=1:nG
    fprintf(1,'%s',SYMBOL{i,1});
    for j=1:nS
        fprintf(1,' %f',OUT(i,j));
    end
    fprintf(1, '\n');
end
```

The loop body is character-for-character identical. The only changes are: `outFile` → `1` (stdout), and the removal of `fopen`/`fclose`. The output format (space-separated, gene symbol first per row) is unchanged.

### 2.3 Why No Other Changes Are Needed

The core numerical algorithm — the `for i = 1:NH` loop, the construction of `i0`, the calls to `quantilenorm` and `CuBlock`, the `log2e` back-conversion, and the assembly of `OUT` — is untouched. The three edits above affect only the I/O boundary of the script, not a single arithmetic operation.

---

## 3. Module Design

### 3.1 `shambhala/io_utils.py` — Unified File I/O

Abstracts the data source (local filesystem vs. S3 object storage) behind a single interface.

**Public functions:**

```python
def read_expression(path: str, sep: str = None) -> pd.DataFrame:
    """
    Read an expression matrix from a local path or S3 URI.

    The file must follow the FL project convention:
        - Rows = samples (row index = sample IDs)
        - Columns = genes (column headers = HGNC gene symbols)
        - First column of the CSV is the sample ID (used as the index)

    Accepts .csv, .tsv, .tsv.gz. Separator is inferred from extension
    if not provided explicitly.

    Parameters
    ----------
    path : str
        Local file path or S3 URI (e.g. 's3://my-bucket/data/input.csv').
    sep : str, optional
        Column separator. Inferred from file extension if None.

    Returns
    -------
    pd.DataFrame
        Shape: (n_samples, n_genes). Index = sample IDs. Columns = gene symbols.
    """
```

```python
def write_expression(df: pd.DataFrame, path: str, sep: str = ',') -> None:
    """
    Write an expression matrix to a local path or S3 URI.

    The DataFrame must follow the FL project convention (rows = samples).
    The index (sample IDs) is written as the first column.

    Parameters
    ----------
    df : pd.DataFrame
        Shape: (n_samples, n_genes).
    path : str
        Destination. Local path or S3 URI.
    sep : str
        Column separator for the output file. Default: ','.
    """
```

```python
def get_s3_client():
    """
    Return a boto3 S3 client using environment credentials or an IAM role.
    Credentials are resolved by the default boto3 credential chain.
    """

def _parse_s3_uri(uri: str) -> tuple[str, str]:
    """
    Parse an S3 URI into (bucket_name, key).

    Example: 's3://my-bucket/path/to/file.csv' → ('my-bucket', 'path/to/file.csv')
    """
```

---

### 3.2 `shambhala/na_handling.py` — NA Strategy

The CuBlock docstring explicitly states: *"The micro array cannot contain NaN values."* Every matrix slice passed to Octave must be NaN-free. This module enforces that guarantee.

**`NAMask` dataclass:**

```python
@dataclass
class NAMask:
    strategy: str                    # 'drop' or 'knn'
    dropped_genes: list[str]         # gene names removed (strategy='drop')
    imputed_positions: pd.DataFrame  # boolean mask of imputed cells (strategy='knn')
    original_gene_order: list[str]   # all gene names in original column order
```

**Public functions:**

```python
def detect_na_genes(df: pd.DataFrame) -> list[str]:
    """
    Return the list of gene names (column names) that contain at least one NaN.

    Parameters
    ----------
    df : pd.DataFrame
        Expression matrix, shape (n_samples, n_genes).

    Returns
    -------
    list[str]
        Gene symbols with one or more NaN values across samples.
    """
```

```python
def apply_na_strategy(
    df: pd.DataFrame,
    strategy: str = 'drop',
    knn_k: int = 5,
    max_na_frac: float = 0.20,
) -> tuple[pd.DataFrame, NAMask]:
    """
    Make an expression matrix NaN-free according to the chosen strategy.

    Strategy 'drop' (default):
        Remove all genes (columns) that contain at least one NaN across any
        sample. The removed genes are recorded in the NAMask and can be
        restored as NaN columns after harmonization via restore_na_genes().

    Strategy 'knn':
        First drop genes where the fraction of NaN values across samples
        exceeds max_na_frac (these are unrecoverable). Then impute remaining
        NaN positions using sklearn.impute.KNNImputer(n_neighbors=knn_k),
        consistent with the KNN imputation used in bench_shared.prepare_dataset_imputed().
        The imputed positions are recorded in the NAMask for logging.
        Imputed values persist in the output (they are not restored to NaN after
        harmonization).

    Parameters
    ----------
    df : pd.DataFrame
        Expression matrix, shape (n_samples, n_genes). May contain NaN.
    strategy : str
        'drop' or 'knn'. Default: 'drop'.
    knn_k : int
        Number of neighbours for KNN imputation. Only used when strategy='knn'.
        Default: 5 (matches bench_shared.prepare_dataset_imputed default).
    max_na_frac : float
        For strategy='knn': genes with a NaN fraction above this threshold are
        dropped before imputation. Default: 0.20 (matches bench_shared default).

    Returns
    -------
    tuple[pd.DataFrame, NAMask]
        (clean_df, mask) where clean_df is NaN-free and mask carries
        metadata needed to reconstruct the original gene set.
    """
```

```python
def restore_na_genes(
    harmonized_df: pd.DataFrame,
    mask: NAMask,
) -> pd.DataFrame:
    """
    Restore dropped genes as NaN columns in the harmonized output.

    Only meaningful when mask.strategy == 'drop'. For strategy 'knn',
    this function returns harmonized_df unchanged (imputed values persist).

    The output preserves the original gene column order from mask.original_gene_order.

    Parameters
    ----------
    harmonized_df : pd.DataFrame
        Harmonized expression matrix without the dropped genes.
    mask : NAMask
        Metadata from apply_na_strategy().

    Returns
    -------
    pd.DataFrame
        Full expression matrix with dropped genes re-inserted as NaN columns,
        in the original gene order.
    """
```

---

### 3.3 `shambhala/octave_bridge.py` — Python → Octave Interface

The core I/O bridge. Formats a genes × samples matrix as TSV, pipes it to Octave via stdin, captures stdout, and parses the result.

**Custom exceptions:**

```python
class OctaveExecutionError(RuntimeError):
    """Raised when Octave exits with a non-zero return code."""

class OctaveOutputError(ValueError):
    """Raised when Octave stdout cannot be parsed or has unexpected shape."""
```

**Public function:**

```python
def run_octave_normalize(
    pool_df: pd.DataFrame,
    nh: int,
    np_: int,
    k: int,
    octave_scripts_dir: str,
    octave_bin: str = 'octave',
    timeout_s: int = 6000,
    random_seed: int = None,
) -> pd.DataFrame:
    """
    Run the full Octave normalization pipeline (quantile normalization + CuBlock)
    on a pool matrix.

    The Octave script Shambhala2_piped.m performs two sequential steps for
    each of the NH input samples (Borisov et al., 2021, doi:10.1093/bioinformatics/btab105;
    Shambhala2.pdf):
        1. Quantile normalization: the input sample is combined with all NP
           P-reference samples and the combined matrix is quantile-normalized
           (Bolstad et al., 2003, doi:10.1093/bioinformatics/19.2.185).
        2. CuBlock normalization: k-means clustering of genes (k clusters) followed
           by cubic polynomial fitting per block, repeated 30 times and averaged.
           Only the normalized values for the NH input samples are retained; the P
           samples serve solely as the quantile-normalization anchor and are discarded.

    This function is the sole interface between Python and Octave. It:
        1. Formats pool_df as a tab-separated string (stdin payload).
        2. Builds the Octave --eval preamble injecting NH, NP, k, and optionally
           a random seed for reproducible k-means initialization.
        3. Calls Octave as a subprocess, piping the TSV via stdin.
        4. Captures and parses Octave's stdout into a DataFrame.
        5. Validates the output shape.

    The pool_df must be in Shambhala2 internal convention:
        - Rows = genes (row index = HGNC gene symbols).
        - Columns = first nh columns are input samples, next np_ columns
          are P reference samples.

    Parameters
    ----------
    pool_df : pd.DataFrame
        Merged input + P reference matrix. Shape: (n_genes, nh + np_).
        Row index = gene symbols (used as the SYMBOL column for Octave).
    nh : int
        Number of input samples in this batch (columns 0..nh-1 of pool_df).
        Within one Octave subprocess, Shambhala2_piped.m loops over each of the
        NH samples individually: each sample is quantile-normalized and
        CuBlock-normalized against the same NP P-reference samples, independently.
        When parallelizing, the full sample set is split into batches; each batch
        runs as a separate subprocess with its own NH. NH=1 is valid.
    np_ : int
        Number of P reference samples (columns nh..nh+np_-1 of pool_df).
    k : int
        Number of k-means probe clusters for CuBlock.
    octave_scripts_dir : str
        Absolute path to the directory containing the .m files.
    octave_bin : str
        Path to the Octave executable. Default: 'octave' (resolved via PATH).
    timeout_s : int
        Subprocess timeout in seconds. Default: 6000 (100 minutes).
    random_seed : int, optional
        If provided, sets Octave's random state via `rand('state', seed)` and
        `randn('state', seed)` in the --eval preamble before sourcing the script.
        Ensures reproducible k-means initialization. Default: None (random).

    Returns
    -------
    pd.DataFrame
        Quantile-normalized and CuBlock-normalized expression values.
        Shape: (n_genes, nh).
        Row index = gene symbols. Columns = input sample names (same as
        the first nh columns of pool_df).

    Raises
    ------
    OctaveExecutionError
        If Octave exits with a non-zero return code. The exception message includes
        the full Octave stderr, a truncated stdout excerpt (first 2000 chars), and
        the exit code, to provide maximum debugging context.
    OctaveOutputError
        If stdout cannot be parsed or the resulting shape does not match
        (n_genes, nh).
    """
```

**Private helpers (with docstrings):**

```python
def _format_pool_as_tsv(pool_df: pd.DataFrame) -> str:
    """
    Serialize pool_df to a tab-separated string in the format expected by
    readExpressionData.m: header row + one row per gene, gene symbol first.

    readExpressionData.m applies log2(x + 1) internally, so values in pool_df
    must be in the raw (non-log) scale. If your data is already log-transformed,
    exponentiate before calling this function.

    Returns
    -------
    str
        Multi-line tab-separated string ready to be piped to Octave stdin.
    """

def _parse_octave_stdout(
    stdout: str,
    expected_genes: list[str],
    expected_samples: list[str],
) -> pd.DataFrame:
    """
    Parse Octave's space-separated stdout into a DataFrame.

    Octave writes one row per gene: '<SYMBOL> <val1> <val2> ... <valN>\n'.
    This function reads that format and returns a properly indexed DataFrame.

    Parameters
    ----------
    stdout : str
        Raw string captured from Octave's stdout.
    expected_genes : list[str]
        Gene symbols in the expected row order (for validation).
    expected_samples : list[str]
        Sample column names to assign to the output DataFrame.

    Returns
    -------
    pd.DataFrame
        Shape: (n_genes, n_samples). Row index = gene symbols.

    Raises
    ------
    OctaveOutputError
        If the parsed shape does not match (len(expected_genes), len(expected_samples)).
    """
```

**Subprocess call pattern:**

```python
seed_preamble = (
    f"rand('state', {random_seed}); randn('state', {random_seed}); "
    if random_seed is not None else ""
)
eval_preamble = (
    f"addpath('{octave_scripts_dir}'); "
    f"{seed_preamble}"
    f"NH={nh}; NP={np_}; k={k}; "
    f"source('Shambhala2_piped.m');"
)
result = subprocess.run(
    [octave_bin, '--no-gui', '--eval', eval_preamble],
    input=tsv_string,
    capture_output=True,
    text=True,
    timeout=timeout_s,
)
if result.returncode != 0:
    raise OctaveExecutionError(
        f"Octave exited with code {result.returncode}.\n"
        f"--- STDERR (full) ---\n{result.stderr}\n"
        f"--- STDOUT (first 2000 chars) ---\n{result.stdout[:2000]}"
    )
```

---

### 3.4 `shambhala/q_rescale.py` — Q-Rescaling Step

Replicates the final rescaling that `Shambhala2.R` performs after Octave returns. This is a direct translation of the R logic into pandas/numpy.

**Algorithm (verbatim from original R):**

```
1. Compute per-gene statistics from Q reference (natural log, not log2):
       rm = mean( log(Q_values + q_pseudocount) )    # per gene, across Q samples
       rs = std(  log(Q_values + q_pseudocount) )    # per gene, across Q samples

2. For each gene g in the Octave output:
       output[g, :] = rm[g] + rs[g] * log(octave_output[g, :] + 1)

3. Exponentiate back:
       output = exp(output)
```

**Public functions:**

```python
def compute_q_statistics(
    q_df: pd.DataFrame,
    q_pseudocount: float = 1e-6,
) -> tuple[pd.Series, pd.Series]:
    """
    Compute per-gene mean and standard deviation of log-expression from the Q
    reference dataset.

    These statistics define the 'universal shape' that every harmonized sample
    will be rescaled to match. Pre-computing them once and passing to rescale()
    avoids recomputing on every parallel batch.

    The original Shambhala2.R applies natural log (not log2) to Q values directly:
        RM = rowMeans(log(Q_matrix))
    This implementation adds q_pseudocount before taking the log:
        rm = mean( log(Q_values + q_pseudocount) )
    making it safe for RNA-seq Q datasets that contain exact zero counts.

    Zero-count handling:
        A pseudocount of 1e-6 (default) is added to all Q values before log.
        This is negligible for non-zero expression values (e.g., log(10 + 1e-6) ≈
        log(10)) but prevents log(0) = -Inf from propagating silently.
        To disable pseudocount and fall back to warn-and-exclude behavior for
        zero-count genes, set q_pseudocount=0: in that case, genes with any zero
        Q value are dropped with a warning, and their statistics are excluded
        from the returned Series.

    Parameters
    ----------
    q_df : pd.DataFrame
        Q reference expression matrix. Shape: (n_samples, n_genes). FL convention.
        Values should be non-negative. Natural log is applied internally after
        adding q_pseudocount.
    q_pseudocount : float
        Small constant added to all Q values before log computation.
        Default: 1e-6. Set to 0 to disable and use warn-and-exclude for zeros instead.

    Returns
    -------
    tuple[pd.Series, pd.Series]
        (rm, rs) — per-gene mean and std of log(Q + q_pseudocount).
        Both indexed by gene symbol. When q_pseudocount=0, genes with any zero
        or negative Q value are excluded.

    Raises
    ------
    ValueError
        If q_df contains negative values (biologically invalid expression levels).
    """

def rescale(
    octave_output: pd.DataFrame,
    rm: pd.Series,
    rs: pd.Series,
) -> pd.DataFrame:
    """
    Apply the Q-based rescaling to the quantile-normalized + CuBlock output.

    For each gene g: output[g] = rm[g] + rs[g] * log(octave_output[g] + 1)
    Then exponentiate: output = exp(output).

    Only genes present in both octave_output and rm/rs are retained
    (inner join on gene symbol). Genes in octave_output but absent from Q
    are dropped with a warning logged.

    Parameters
    ----------
    octave_output : pd.DataFrame
        Quantile-normalized and CuBlock-normalized output from Octave.
        Shape: (n_genes, n_samples). Rows = genes, columns = samples.
    rm : pd.Series
        Per-gene mean of log(Q + pseudocount), indexed by gene symbol.
    rs : pd.Series
        Per-gene std of log(Q + pseudocount), indexed by gene symbol.

    Returns
    -------
    pd.DataFrame
        Fully harmonized expression matrix. Same shape as octave_output
        (after inner join on gene set). Rows = genes, columns = samples.
    """
```

---

### 3.5 `shambhala/parallel.py` — Parallel Sample Dispatcher

Splits samples into batches and dispatches each batch to an isolated worker process.

**Why `ProcessPoolExecutor`, not `ThreadPoolExecutor`:**
Each worker calls `subprocess.run` (Octave). Multiple threads sharing the same Python process can have file descriptor conflicts on `stdin`/`stdout` pipes. Process-level isolation gives each worker completely independent file descriptors, environment, and memory space.

**Worker function (module-level, required for pickling):**

```python
def _harmonize_one_batch(
    batch_sample_names: list[str],
    input_df: pd.DataFrame,
    p_df: pd.DataFrame,
    k: int,
    octave_scripts_dir: str,
    octave_bin: str,
    timeout_s: int,
    random_seed: int,
) -> pd.DataFrame:
    """
    Harmonize one batch of samples through the full Octave normalization pipeline.

    This function is executed in a worker process. It must be module-level
    (not a lambda or closure) to be picklable by multiprocessing.

    Steps:
        1. Extract the batch columns from input_df.
        2. Concatenate with the full P reference: pool = [batch | P].
        3. Format pool as TSV and call run_octave_normalize().
        4. Return the quantile-normalized + CuBlock-normalized batch (genes × batch_samples).

    Parameters
    ----------
    batch_sample_names : list[str]
        Sample IDs (column names in input_df) assigned to this batch.
    input_df : pd.DataFrame
        Full input matrix. Shape: (n_genes, n_samples). Rows = genes.
    p_df : pd.DataFrame
        P reference matrix. Shape: (n_genes, n_p_samples). Rows = genes.
    k : int
        k-means clusters for CuBlock.
    octave_scripts_dir : str
        Absolute path to the octave/ directory.
    octave_bin : str
        Path to the Octave executable.
    timeout_s : int
        Per-batch Octave subprocess timeout in seconds.
    random_seed : int
        Passed through to run_octave_normalize for reproducible k-means.

    Returns
    -------
    pd.DataFrame
        Quantile-normalized + CuBlock-normalized batch.
        Shape: (n_genes, len(batch_sample_names)).
    """
```

**Public dispatcher function:**

```python
def harmonize_parallel(
    input_df: pd.DataFrame,
    p_df: pd.DataFrame,
    rm: pd.Series,
    rs: pd.Series,
    k: int = 5,
    n_workers: int = 5,
    octave_scripts_dir: str = None,
    octave_bin: str = 'octave',
    timeout_s: int = 6000,
    random_seed: int = None,
) -> pd.DataFrame:
    """
    Harmonize all samples in parallel and apply Q-rescaling to the assembled result.

    Samples are split into n_workers batches of approximately equal size using
    numpy.array_split. Each batch is dispatched to a separate worker process
    that runs an independent Octave subprocess. Results are concatenated and
    Q-rescaling is applied once to the full assembled matrix.

    Parameters
    ----------
    input_df : pd.DataFrame
        NaN-free expression matrix in internal Shambhala convention.
        Shape: (n_genes, n_samples). Rows = genes, columns = samples.
    p_df : pd.DataFrame
        NaN-free P reference. Shape: (n_genes, n_p_samples).
    rm : pd.Series
        Per-gene Q mean (from q_rescale.compute_q_statistics).
    rs : pd.Series
        Per-gene Q std.
    k : int
        k-means clusters for CuBlock. Default: 5.
    n_workers : int
        Number of parallel worker processes. Default: 5. Maximum: 10.
        If n_workers > n_samples, it is clamped to n_samples.
    octave_scripts_dir : str
        Absolute path to the directory containing the .m files.
        If None, defaults to the `octave/` directory adjacent to this module.
    octave_bin : str
        Path to the Octave executable. Default: 'octave'.
    timeout_s : int
        Per-batch Octave subprocess timeout in seconds. Default: 6000.
    random_seed : int, optional
        If provided, passed to every worker for reproducible k-means.
        Default: None (random initialization).

    Returns
    -------
    pd.DataFrame
        Fully harmonized expression matrix. Shape: (n_genes, n_samples).
        Rows = genes, columns = samples (same order as input).

    Raises
    ------
    ValueError
        If n_workers is outside the range [1, 10].
    RuntimeError
        If any worker process raises an exception (original exception chained).
    """
```

**Parallelization logic:**

```
samples = list(input_df.columns)                      # n_samples column names
batches = numpy.array_split(samples, n_workers)       # split into n_workers groups

with ProcessPoolExecutor(max_workers=n_workers) as executor:
    futures = {
        executor.submit(_harmonize_one_batch, batch, input_df, p_df, ...): batch
        for batch in batches
    }
    results = []
    for future in as_completed(futures):
        results.append(future.result())               # raises if worker failed

assembled = pd.concat(results, axis=1)[samples]       # restore original column order
harmonized = q_rescale.rescale(assembled, rm, rs)
return harmonized
```

---

### 3.6 `run_shambhala.py` — CLI Entry Point

This is the single file a user or pipeline script ever needs to call.

#### Command-Line Interface

```
usage: run_shambhala.py [-h]
                        --input PATH
                        --P PATH
                        --Q PATH
                        --output PATH
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

Harmonize a gene expression matrix using the Shambhala2 algorithm
(quantile normalization + CuBlock + Q-rescaling, one sample at a time).

required arguments:
  --input PATH          Input expression matrix.
                        Local path or S3 URI (e.g. s3://bucket/input.csv).
                        Convention: rows=samples, columns=genes.
                        Formats: .csv, .tsv, .tsv.gz (auto-detected).
  --P PATH              Calibration reference dataset P.
                        Local path or S3 URI. Same format as --input.
                        P is used as the quantile-normalization anchor: each input
                        sample is jointly quantile-normalized with all P samples.
  --Q PATH              Definitive reference dataset Q.
                        Local path or S3 URI. Same format as --input.
                        Q defines the universal expression shape: per-gene mean
                        and std of log(Q + pseudocount) are used to rescale output.
  --output PATH         Destination for the harmonized output matrix.
                        Local path or S3 URI. Written in CSV format.

optional arguments:
  --k INT               Number of k-means gene clusters for the CuBlock algorithm.
                        Default: 5. Do not change without benchmarking — this
                        parameter governs the granularity of the piecewise-cubic
                        fitting inside Octave.
  --n-workers INT       Number of parallel worker processes. Each worker runs an
                        independent Octave subprocess on a separate batch of samples.
                        Default: 5. Range: [1, 10].
  --na-strategy {drop,knn}
                        Strategy for handling NaN values in the input matrix.
                        CuBlock requires NaN-free data; this strategy is applied
                        before passing data to Octave.
                          drop (default): Remove genes with any NaN across samples.
                            Dropped genes are restored as NaN columns in the output.
                          knn: Impute NaN values using KNN (sklearn KNNImputer).
                            Genes where NaN fraction exceeds --max-na-frac are
                            dropped first. Imputed values persist in the output.
                        Default: drop.
  --knn-k INT           Number of neighbours for KNN imputation.
                        Only used when --na-strategy knn. Default: 5.
  --max-na-frac FLOAT   Maximum allowed NaN fraction per gene before dropping.
                        Only used when --na-strategy knn. Default: 0.20.
  --q-pseudocount FLOAT Small constant added to all Q values before taking log,
                        to prevent log(0) = -Inf for zero-count RNA-seq Q datasets.
                        Default: 1e-6. Set to 0 to disable (genes with zero counts
                        in Q are then excluded with a warning instead).
  --octave-bin PATH     Path to the Octave executable.
                        Default: 'octave' (resolved via PATH).
                        For K8s pods: /opt/conda/bin/octave.
  --timeout-s INT       Per-batch Octave subprocess timeout in seconds.
                        Default: 6000. Increase for very large batches or slow nodes.
  --random-seed INT     Integer seed for Octave's random state (rand and randn).
                        Ensures reproducible k-means initialization in CuBlock.
                        Default: not set (random initialization).
  --log-level {DEBUG,INFO,WARNING}
                        Logging verbosity. Default: INFO.
```

#### Full Execution Flow

```
Step 1  — Parse and validate arguments.
          Enforce n_workers ∈ [1, 10].

Step 2  — Read input:
          input_df = io_utils.read_expression(args.input)      # (n_samples, n_genes)
          p_df     = io_utils.read_expression(args.P)          # (n_p_samples, n_genes)
          q_df     = io_utils.read_expression(args.Q)          # (n_q_samples, n_genes)
          Log: shape of each matrix.

Step 3  — Transpose to internal convention (genes × samples):
          input_T = input_df.T   # (n_genes, n_samples)
          p_T     = p_df.T       # (n_genes, n_p_samples)
          q_T     = q_df.T       # (n_genes, n_q_samples)

Step 4  — Find common genes across input, P, and Q (inner join).
          Log the count of genes retained and dropped per dataset.

Step 5  — Apply NA strategy to the input matrix only:
          (clean_T, na_mask) = na_handling.apply_na_strategy(
              input_T, args.na_strategy, knn_k=args.knn_k,
              max_na_frac=args.max_na_frac,
          )
          Log: n_genes affected, strategy used.

Step 6  — Pre-compute Q statistics (done once, reused by all workers):
          (rm, rs) = q_rescale.compute_q_statistics(
              q_T, q_pseudocount=args.q_pseudocount,
          )
          Log: n_genes with valid Q statistics; n_genes excluded if any.

Step 7  — Dispatch parallel harmonization:
          harmonized_T = parallel.harmonize_parallel(
              input_df=clean_T, p_df=p_T, rm=rm, rs=rs,
              k=args.k, n_workers=args.n_workers,
              octave_scripts_dir=<abs path to octave/>,
              octave_bin=args.octave_bin,
              timeout_s=args.timeout_s,
              random_seed=args.random_seed,
          )

Step 8  — Restore NA genes (if strategy='drop'):
          full_harmonized_T = na_handling.restore_na_genes(harmonized_T, na_mask)

Step 9  — Transpose back to FL convention (samples × genes):
          result_df = full_harmonized_T.T   # (n_samples, n_genes)

Step 10 — Write output:
          io_utils.write_expression(result_df, args.output)

Step 11 — Log summary:
          - Total samples processed
          - Total genes in output
          - NA strategy applied and genes affected
          - n_workers used
          - Wall-clock time elapsed
```

#### Logging convention

All steps emit structured log messages at `INFO` level. Worker-level errors include the batch sample list in the message. Overall progress uses a running counter:
```
[INFO] Step 2/11 — Reading input: s3://bucket/input.csv (7237 samples × 3447 genes)
[INFO] Step 7/11 — Dispatching 7237 samples across 5 workers (1448 samples/worker)
[INFO]   Worker 1/5 complete: 1448 samples in 142.3 s
...
```

---

## 4. Data Convention Details

### 4.1 External convention (FL project, all file I/O)

```
rows    = samples  (index = sample IDs, e.g. GSM123456)
columns = genes    (headers = HGNC symbols, e.g. TP53, CD19)
```

This is what `read_expression` returns and `write_expression` expects. All fixture test files are in this format.

### 4.2 Internal convention (inside Octave bridge and Q-rescaling)

```
rows    = genes    (index = HGNC symbols)
columns = samples  (headers = sample IDs)
```

This matches what `readExpressionData.m` expects (first column = SYMBOL, remaining columns = samples). The transpose happens in `run_shambhala.py` at Steps 3 and 9.

### 4.3 Pool matrix for Octave (genes × [input_batch + P])

```
rows    = genes (genes present in intersection of input, P, Q)
columns = [ sample_batch_1 | sample_batch_2 | ... | P_1 | P_2 | ... ]
            ← nh columns →                         ← np_ columns →
```

The SYMBOL column (gene symbols) is prepended from the DataFrame row index when formatting the TSV for stdin.

---

## 5. Testing Plan

### 5.1 Fixture Preparation

`tests/generate_fixtures.py` generates all fixture files from the existing `Shambhala2/` data. Run once before executing the test suite:

```bash
cd Shambhala_containerized/
python tests/generate_fixtures.py
```

Files produced:
- `input_10samples.csv` — transpose `Shambhala2/Input.csv` (genes → columns, samples → rows)
- `expected_output.csv` — transpose `Shambhala2/Output.csv` the same way
- `P0_small.csv` — transpose `Shambhala2/P0.csv` (Affymetrix GPL570 calibration reference)
- `Q0_small.csv` — transpose `Shambhala2/Q0.csv` (GTEx definitive reference)
- `input_with_nas.csv` — copy of `input_10samples.csv` with 5% of values replaced by NaN (`numpy.random.seed(42)`)

### 5.2 `test_integration.py` — End-to-End Tests

The full test matrix below is also documented in `README.md` for pre-run review.

| Test | Description | Pass Criterion |
|---|---|---|
| `test_single_worker_matches_reference` | Run with `n_workers=1` on 10-sample fixture | Output matches `expected_output.csv` within `atol=1e-4` |
| `test_parallel_matches_single` | Run with `n_workers=3` and `n_workers=5` | Both match the `n_workers=1` result within `atol=0.05` on per-gene row means (stochastic k-means) |
| `test_local_io_roundtrip` | Write to temp file; read back | Shape and values identical |
| `test_gene_intersection_is_applied` | P has 50 fewer genes than input | Output gene count = intersection size |
| `test_na_drop_full_pipeline` | Run on `input_with_nas.csv` with `--na-strategy drop` | Output has NaN exactly where NA genes were; non-NA genes match expected values |
| `test_na_knn_full_pipeline` | Run on `input_with_nas.csv` with `--na-strategy knn` | Output has no NaN; shape equals full input shape |

### 5.3 `test_na_handling.py` — Unit Tests

| Test | Description |
|---|---|
| `test_detect_na_genes_correct` | Input with 3 known NaN genes → returns exactly those 3 |
| `test_drop_removes_na_genes` | After `apply_na_strategy(..., 'drop')`, returned df has zero NaN |
| `test_knn_fills_nas` | After `apply_na_strategy(..., 'knn')`, returned df has zero NaN |
| `test_restore_preserves_original_order` | After `restore_na_genes`, column order matches original |
| `test_restore_inserts_nan_correct_columns` | After `restore_na_genes`, NaN is in exactly the dropped gene columns |
| `test_all_na_gene_handled_gracefully` | Gene with 100% NaN under 'knn' (exceeds max_na_frac) → dropped, not imputed; no exception |

### 5.4 `test_octave_bridge.py` — Unit Tests

| Test | Description |
|---|---|
| `test_bridge_returns_correct_shape` | 5-sample batch + P → output shape is (n_genes, 5) |
| `test_bridge_symbol_index_correct` | Row index of output matches gene symbols from input |
| `test_bridge_raises_on_nonzero_exit` | Pass invalid data → `OctaveExecutionError` raised |
| `test_bridge_error_message_contains_stderr` | Trigger Octave error → exception message includes full Octave stderr |
| `test_bridge_raises_on_unparseable_stdout` | Patch subprocess to return garbage stdout → `OctaveOutputError` raised |
| `test_tsv_formatting_round_trip` | `_format_pool_as_tsv` then `pd.read_csv` → values match original |
| `test_random_seed_injected_in_preamble` | With `random_seed=42` → `--eval` string contains `rand('state', 42)` |

### 5.5 `test_q_rescale.py` — Unit Tests

| Test | Description |
|---|---|
| `test_rescale_known_values` | Toy input with hand-computed expected output → matches within 1e-10 |
| `test_compute_q_statistics_shape` | Returns two Series with index = gene symbols |
| `test_pseudocount_prevents_log_zero` | Q with a zero in one gene + default pseudocount → no exclusion, no warning |
| `test_pseudocount_zero_excludes_with_warning` | Q with a zero in one gene + `q_pseudocount=0` → gene excluded, warning logged |
| `test_rescale_raises_on_negative_q` | Q with a negative value → `ValueError` from `compute_q_statistics` |
| `test_rescale_gene_intersection` | Input has 5 genes, Q has 4 → output has 4 genes |

### 5.6 `test_io_utils.py` — Unit Tests

| Test | Description |
|---|---|
| `test_read_local_csv` | Read fixture CSV → correct shape and dtypes |
| `test_read_local_tsv_gz` | Read gzipped TSV → same result as uncompressed |
| `test_write_and_read_roundtrip` | Write df to temp file; read back → identical |
| `test_parse_s3_uri` | `'s3://b/k/f.csv'` → `('b', 'k/f.csv')` |
| `test_s3_read_mocked` | Mock `boto3.client.get_object`; assert `read_expression` calls it with correct bucket/key |

---

## 6. Dependencies

```
# Runtime — pure Python, no R, no rpy2
pandas >= 1.5
numpy  >= 1.23, < 2.0          # NumPy < 2.0 per CLAUDE.md requirement
scikit-learn >= 1.0             # KNNImputer for --na-strategy knn
boto3  >= 1.26                  # S3 I/O (optional; only needed for S3 paths)

# System
octave                          # Any version supporting /dev/stdin and fprintf(1,...)
                                # Tested with Conda Octave at /opt/conda/bin/octave

# Dev / testing only
pytest
pytest-timeout                  # Prevents hanging Octave subprocesses in CI
```

No `rpy2`. No R. The R logic from `Shambhala2.R` is replaced by a few lines of pandas/numpy in the `q_rescale` module. This simplifies the dependency stack substantially and is the principal enabler of clean multi-process execution.

---

## 7. Resolved Design Decisions

All questions were reviewed and approved by Daniil.

| Decision | Resolution |
|---|---|
| R elimination | Approved. `Shambhala_containerized/` is fully independent; `Shambhala2/` and `bench_shared.py` are untouched. |
| Input format for P and Q | FL convention throughout (rows = samples, columns = genes). The script transposes internally. |
| Gene name column | Columns = gene symbols; first CSV column = sample ID (used as index). No dedicated `SYMBOL` column in input files. |
| NA strategy default | `drop` is the default. `knn` is the alternative, activated via `--na-strategy knn`. Median imputation is not used. |
| KNN imputation implementation | `sklearn.impute.KNNImputer`, consistent with `bench_shared.prepare_dataset_imputed`. Default `knn_k=5`, `max_na_frac=0.20`. |
| Q pseudocount | `q_pseudocount=1e-6` (default): added to all Q values before log to handle RNA-seq zero counts safely. Set `--q-pseudocount 0` to disable and fall back to warn-and-exclude. |
| Octave binary default | `'octave'` (resolved via PATH). K8s invocations pass `--octave-bin /opt/conda/bin/octave` explicitly. |
| Parallel reproducibility | `--random-seed` CLI argument seeds Octave's `rand` and `randn` state via the `--eval` preamble. |
| Timeout default | `--timeout-s` defaults to 6000 seconds. |

---

## 8. What Is Explicitly Out of Scope

- No changes to `Shambhala2/` (original directory is untouched).
- No changes to `harmonization-scripts/bench_shared.py` or `run_norm_job.py`.
- No Dockerfile or K8s manifest (the existing `harmonization-scripts/k8s/CLAUDE.md` covers infrastructure; environment setup is out of scope for this plan).
- No modification of P/Q reference datasets (use the existing `Shambhala2/P0.csv` and `Q0.csv`).
- No change to the CuBlock algorithm, quantile normalization, or k-means initialization strategy (these are inherited verbatim from `Shambhala2/`).

---

## 9. Implementation To-Do List

### Phase 0 — Repository Scaffolding

- [x] Create `Shambhala_containerized/` directory
- [x] Create `shambhala/__init__.py` (empty; marks the package)
- [x] Create `octave/`, `tests/`, `tests/fixtures/` directories
- [x] Write `pyproject.toml` — package name `shambhala-containerized`, version `0.1.0`, `python_requires = ">=3.11"`, authors, license
- [x] Write `requirements.txt` — pin `pandas`, `numpy<2.0`, `scikit-learn`, `boto3`, `pytest`, `pytest-timeout`
- [x] Write `setup.py` — minimal shim calling `setuptools.setup()` for `pip install -e .` compatibility

### Phase 1 — Octave Scripts

- [x] Copy `Shambhala2/CuBlock.m` → `octave/CuBlock.m` verbatim (verified byte-for-byte identical)
- [x] Copy `Shambhala2/quantilenorm.m` → `octave/quantilenorm.m` verbatim
- [x] Copy `Shambhala2/kmeans.m` → `octave/kmeans.m` verbatim
- [x] Copy `Shambhala2/readExpressionData.m` → `octave/readExpressionData.m` verbatim
- [x] Create `octave/Shambhala2_piped.m` by applying exactly 3 edits to `Shambhala2/Shambhala2.m`:
  - [x] Edit 1: delete the 5-line `args.txt` reading block
  - [x] Edit 2: replace `readExpressionData("P_prim.txt",'log2')` with `readExpressionData('/dev/stdin','log2')`
  - [x] Edit 3: replace the `fopen/fprintf-to-file/fclose` output block with `fprintf(1,...)` loop
- [x] Run `diff Shambhala2/Shambhala2.m octave/Shambhala2_piped.m` and confirm exactly 3 change sites

### Phase 2 — Test Fixture Generation

- [x] Write `tests/generate_fixtures.py`:
  - [x] Read `Shambhala2/Input.csv` and transpose → `fixtures/input_10samples.csv`
  - [x] Read `Shambhala2/Output.csv` and transpose → `fixtures/expected_output.csv`
  - [x] Read `Shambhala2/P0.csv` and transpose → `fixtures/P0_small.csv`
  - [x] Read `Shambhala2/Q0.csv` and transpose → `fixtures/Q0_small.csv`
  - [x] Copy `input_10samples.csv`, inject 5% NaN with `numpy.random.seed(42)` → `fixtures/input_with_nas.csv`
- [x] Run `python tests/generate_fixtures.py` and verify all 5 files are created with correct shapes

### Phase 3a — `shambhala/io_utils.py`

- [x] Implement `_parse_s3_uri(uri)` → `(bucket, key)`
- [x] Implement `get_s3_client()` using boto3 default credential chain
- [x] Implement `read_expression(path, sep=None)`:
  - [x] Detect `s3://` prefix; download via `get_object` into `BytesIO` buffer
  - [x] Infer separator from file extension (`.csv` → `,`, `.tsv` / `.tsv.gz` → `\t`)
  - [x] Handle `.gz` decompression automatically via `pd.read_csv(..., compression='gzip')`
  - [x] Set first column as index (sample IDs)
- [x] Implement `write_expression(df, path, sep=',')`:
  - [x] Write index as first column
  - [x] S3 path: serialize to `StringIO`, call `put_object`

### Phase 3b — `shambhala/na_handling.py`

- [x] Define `NAMask` dataclass with fields: `strategy`, `dropped_genes`, `imputed_positions`, `original_gene_order`
- [x] Implement `detect_na_genes(df)` → `list[str]`
- [x] Implement `apply_na_strategy` — `'drop'` branch:
  - [x] Identify all genes with ≥1 NaN
  - [x] Drop them; record in `NAMask.dropped_genes`
  - [x] Store original column order in `NAMask.original_gene_order`
- [x] Implement `apply_na_strategy` — `'knn'` branch:
  - [x] Drop genes where `na_frac > max_na_frac`; record in `NAMask.dropped_genes`
  - [x] Fit `KNNImputer(n_neighbors=knn_k)` and transform
  - [x] Record imputed positions (boolean mask) in `NAMask.imputed_positions`
- [x] Implement `restore_na_genes(harmonized_df, mask)`:
  - [x] Return unchanged if `mask.strategy == 'knn'`
  - [x] Re-insert dropped genes as NaN columns in original order

### Phase 3c — `shambhala/q_rescale.py`

- [x] Implement `compute_q_statistics(q_df, q_pseudocount=1e-6)`:
  - [x] Raise `ValueError` for any negative Q values
  - [x] When `q_pseudocount > 0`: compute `np.log(q_T + q_pseudocount)` for all genes
  - [x] When `q_pseudocount == 0`: detect zero-count genes, log warning, exclude them; compute `np.log(q_T)` on remaining
  - [x] Return `(rm, rs)` — per-gene mean and std, indexed by gene symbol
- [x] Implement `rescale(octave_output, rm, rs)`:
  - [x] Inner join on gene symbol between `octave_output` index and `rm`/`rs` index
  - [x] Warn and log count of genes dropped by the inner join
  - [x] Apply: `result = np.exp(rm + rs * np.log(octave_output + 1))`
  - [x] Return DataFrame with gene rows and sample columns

### Phase 3d — `shambhala/octave_bridge.py`

- [x] Define `OctaveExecutionError(RuntimeError)` and `OctaveOutputError(ValueError)`
- [x] Implement `_format_pool_as_tsv(pool_df)`:
  - [x] Header row: `SYMBOL\tsample1\tsample2\t...`
  - [x] One row per gene: `GENESYM\tval1\tval2\t...`
  - [x] Return as a single string (no temporary files)
- [x] Implement `_parse_octave_stdout(stdout, expected_genes, expected_samples)`:
  - [x] Parse space-separated rows: first token = gene symbol, rest = float values
  - [x] Validate shape against `(len(expected_genes), len(expected_samples))`
  - [x] Raise `OctaveOutputError` on mismatch
- [x] Implement `run_octave_normalize(pool_df, nh, np_, k, octave_scripts_dir, octave_bin, timeout_s, random_seed)`:
  - [x] Build `eval_preamble` with optional seed injection
  - [x] Call `subprocess.run` with `input=tsv_string`, `capture_output=True`
  - [x] On non-zero exit code: raise `OctaveExecutionError` with full stderr + stdout excerpt
  - [x] Parse stdout via `_parse_octave_stdout`; return DataFrame

### Phase 3e — `shambhala/parallel.py`

- [x] Implement module-level `_harmonize_one_batch(batch_sample_names, input_df, p_df, k, octave_scripts_dir, octave_bin, timeout_s, random_seed)`:
  - [x] Extract batch columns; concatenate with p_df: `pool = pd.concat([batch, p_df], axis=1)`
  - [x] Call `run_octave_normalize(pool, nh=len(batch), np_=p_df.shape[1], ...)`
  - [x] Return normalized batch DataFrame
- [x] Implement `harmonize_parallel(...)`:
  - [x] Validate `n_workers ∈ [1, 10]`; clamp if `n_workers > n_samples`
  - [x] Resolve `octave_scripts_dir` default to `octave/` adjacent to module
  - [x] Split samples with `numpy.array_split`
  - [x] Dispatch with `ProcessPoolExecutor`; collect results with `as_completed`
  - [x] Re-raise worker exceptions with original traceback chained
  - [x] Concatenate results in original sample order
  - [x] Apply `q_rescale.rescale(assembled, rm, rs)`

### Phase 4 — CLI (`run_shambhala.py`)

- [x] Set up `argparse.ArgumentParser` with all 13 arguments (required + optional) per §3.6
- [x] Add validation: `n_workers` must be in [1, 10]
- [x] Configure `logging` at the requested level; use step-counter format
- [x] Implement all 11 steps from §3.6:
  - [x] Steps 1–4: argument parsing, reading, transposing, gene intersection
  - [x] Step 5: NA strategy application
  - [x] Step 6: Q statistics with `q_pseudocount`
  - [x] Step 7: parallel harmonization dispatch
  - [x] Steps 8–11: restore NAs, transpose back, write, log summary
- [x] Resolve absolute path to `octave/` directory relative to this script's location

### Phase 5 — Test Suite

- [x] Write `tests/test_io_utils.py` — 5 tests (§5.6); all 5 pass
- [x] Write `tests/test_na_handling.py` — 6 tests (§5.3); all 6 pass
- [x] Write `tests/test_q_rescale.py` — 6 tests (§5.5); all 6 pass
- [x] Write `tests/test_octave_bridge.py` — 7 tests (§5.4); 5 pass, 2 skip (require Octave)
- [x] Write `tests/test_integration.py` — 6 tests (§5.2); 0 pass / 6 skip (require Octave)
- [x] Run full suite: `pytest tests/ -v` — 22 passed, 8 skipped (Octave not installed on this machine; skipped tests have `@octave_required` decorator and will run when Octave is present)

### Phase 6 — Documentation

- [x] Write `README.md`:
  - [x] Algorithm section: three-step Shambhala2 process with citations
  - [x] Installation section: `pip install -e .`; system prerequisites (Octave, Python 3.11)
  - [x] Full CLI reference section: copy the §3.6 usage block with two worked examples
  - [x] Testing section: reproduce the test tables from §5.2–5.6 so Daniil can review coverage
- [x] Update `shambhala_adoption/CLAUDE.md`: mark `Shambhala_containerized/` as implemented

### Phase 7 — Smoke Test on Real FL Data

- [ ] Run on 100-sample subset to validate against previous R-pipeline output:
  ```bash
  cd Shambhala_containerized/
  python run_shambhala.py \
      --input ../Shambhala2/large_tables/comb_exp_log_dedup_strict_subset.csv \
      --P ../Shambhala2/P0.csv \
      --Q ../Shambhala2/Q0.csv \
      --output /tmp/shambhala_containerized_subset.csv \
      --n-workers 5 --random-seed 42 --log-level DEBUG
  ```
  **Note:** Requires Octave. The fixture input files (`P0.csv`, `Q0.csv`) are in Shambhala2-internal convention (genes × samples); use `--P ../Shambhala2/P0.csv` directly — `run_shambhala.py` reads these with `io_utils.read_expression` which expects FL convention (samples × genes). Use the transposed copies in `tests/fixtures/P0_small.csv` and `tests/fixtures/Q0_small.csv` for the smoke test.
- [ ] Compare `/tmp/shambhala_containerized_subset.csv` against `Shambhala2/large_tables/shambhala_harmonized_strict_subset.csv` (ground truth from original R pipeline) — per-gene row mean should agree within `atol=0.1`
- [ ] Document observed differences (expected: minor floating-point divergence from random k-means initialization) in a comment in `CLAUDE.md`
