"""
Python → Octave subprocess interface for the Shambhala2 normalization pipeline.

This module is the sole interface between Python and Octave. It formats data
as TSV, pipes it to Octave via stdin, and parses the space-separated stdout.
No intermediate files are written.
"""

import re
import subprocess
import tempfile
import threading
import time
from io import StringIO

import numpy as np
import pandas as pd

from shambhala.progress_display import ProgressEvent, _QueueLike


class OctaveExecutionError(RuntimeError):
    """Raised when Octave exits with a non-zero return code."""


class OctaveOutputError(ValueError):
    """Raised when Octave stdout cannot be parsed or has unexpected shape."""


def _format_pool_as_tsv(pool_df: pd.DataFrame) -> str:
    """
    Serialize pool_df to a tab-separated string in the format expected by
    readExpressionData.m: header row + one row per gene, gene symbol first.

    readExpressionData.m applies log2(x + 1) internally, so values in pool_df
    must be in the raw (non-log) scale. If your data is already log-transformed,
    exponentiate before calling this function.

    Parameters
    ----------
    pool_df : pd.DataFrame
        Genes × samples matrix. Row index = gene symbols.
        Columns = sample IDs.

    Returns
    -------
    str
        Multi-line tab-separated string ready to be piped to Octave stdin.
    """
    buf = StringIO()
    # Header: SYMBOL <tab> sample1 <tab> sample2 ...
    buf.write("SYMBOL\t" + "\t".join(str(c) for c in pool_df.columns) + "\n")
    for gene_sym, row in pool_df.iterrows():
        buf.write(str(gene_sym) + "\t" + "\t".join(str(v) for v in row.values) + "\n")
    return buf.getvalue()


def _parse_octave_stdout(
    stdout: str,
    expected_genes: list[str],
    expected_samples: list[str],
) -> pd.DataFrame:
    """
    Parse Octave's space-separated stdout into a DataFrame.

    Octave writes one row per gene: '<SYMBOL> <val1> <val2> ... <valN>\\n'.
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
    rows = []
    gene_names = []
    for line in stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        tokens = line.split()
        if len(tokens) < 2:
            continue
        gene_sym = tokens[0]
        try:
            values = [float(t) for t in tokens[1:]]
        except ValueError:
            # Non-data lines (e.g. Octave disp() messages) — skip
            continue
        if len(values) == len(expected_samples):
            gene_names.append(gene_sym)
            rows.append(values)

    if len(rows) != len(expected_genes) or len(gene_names) != len(expected_genes):
        raise OctaveOutputError(
            f"Parsed {len(rows)} gene rows from Octave stdout, "
            f"expected {len(expected_genes)}. "
            f"First 500 chars of stdout: {stdout[:500]!r}"
        )

    return pd.DataFrame(rows, index=gene_names, columns=expected_samples)


def run_octave_normalize(
    pool_df: pd.DataFrame,
    nh: int,
    np_: int,
    k: int,
    octave_scripts_dir: str,
    octave_bin: str = "octave",
    timeout_s: int = 6000,
    random_seed: int | None = None,
    fixed_clusters: np.ndarray | None = None,
    skip_qn: bool = False,
) -> pd.DataFrame:
    """
    Run the full Octave normalization pipeline (quantile normalization + CuBlock)
    on a pool matrix.

    The Octave script Shambhala2_piped.m performs two sequential steps for
    each of the NH input samples (Borisov et al., 2021, doi:10.1093/bioinformatics/btab105):
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
        NH=1 is valid. When parallelizing, the full sample set is split into
        batches; each batch runs as a separate subprocess with its own NH.
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
    input_sample_names = list(pool_df.columns[:nh])
    gene_symbols = list(pool_df.index)

    tsv_string = _format_pool_as_tsv(pool_df)

    seed_preamble = (
        f"rand('state', {random_seed}); randn('state', {random_seed}); "
        if random_seed is not None
        else ""
    )

    _clusters_tmpfile_sync = None
    clusters_preamble_sync = ""
    if skip_qn:
        script_name_sync = "Shambhala2_piped_preqn.m"
    else:
        script_name_sync = "Shambhala2_piped.m"
    if fixed_clusters is not None:
        _clusters_tmpfile_sync = tempfile.NamedTemporaryFile(
            mode="w", suffix=".txt", delete=False
        )
        for label in fixed_clusters:
            _clusters_tmpfile_sync.write(f"{int(label)}\n")
        _clusters_tmpfile_sync.flush()
        _clusters_tmpfile_sync.close()
        clusters_preamble_sync = (
            f"FIXED_CLUSTERS = int32(load('{_clusters_tmpfile_sync.name}')); "
        )
        script_name_sync = "Shambhala2_piped_fixed.m"

    eval_preamble = (
        f"addpath('{octave_scripts_dir}'); "
        f"{seed_preamble}"
        f"NH={nh}; NP={np_}; k={k}; "
        f"{clusters_preamble_sync}"
        f"source('{octave_scripts_dir}/{script_name_sync}');"
    )

    try:
        result = subprocess.run(
            [octave_bin, "--no-gui", "--eval", eval_preamble],
            input=tsv_string,
            capture_output=True,
            text=True,
            timeout=timeout_s,
        )
    except FileNotFoundError:
        raise OctaveExecutionError(
            f"Octave executable not found: {octave_bin!r}. "
            "Install Octave or pass the correct path via --octave-bin."
        )

    if _clusters_tmpfile_sync is not None:
        import os as _os  # noqa: PLC0415
        try:
            _os.unlink(_clusters_tmpfile_sync.name)
        except OSError:
            pass

    if result.returncode != 0:
        raise OctaveExecutionError(
            f"Octave exited with code {result.returncode}.\n"
            f"--- STDERR (full) ---\n{result.stderr}\n"
            f"--- STDOUT (first 2000 chars) ---\n{result.stdout[:2000]}"
        )

    return _parse_octave_stdout(result.stdout, gene_symbols, input_sample_names)


def precompute_clusters_from_p(
    p_T: pd.DataFrame,
    k: int = 5,
    random_seed: int | None = None,
) -> np.ndarray:
    """
    Precompute stable gene cluster assignments from the P reference data.

    Runs sklearn KMeans on log2(P+1) with n_init=30 to get a stable partition
    for use with CuBlock_fixed.m. Clustering is on genes (rows), using P sample
    expression profiles (columns) as features.

    Parameters
    ----------
    p_T : pd.DataFrame
        P reference matrix in internal convention: rows=genes, columns=P samples.
    k : int
        Number of clusters. Must match the k used in the main pipeline. Default 5.
    random_seed : int, optional
        Seed for reproducible KMeans initialization.

    Returns
    -------
    np.ndarray
        1-indexed cluster labels (1..k), shape (G,), dtype int32. Ready for
        injection into Octave as FIXED_CLUSTERS.
    """
    from sklearn.cluster import KMeans  # noqa: PLC0415

    p_log2 = np.log2(p_T.values + 1.0)  # (G, NP) — match Octave's log2(x+1)
    km = KMeans(
        n_clusters=k,
        n_init=30,
        max_iter=1000,
        random_state=random_seed,
    )
    labels_0indexed = km.fit_predict(p_log2)  # cluster genes (rows); features = P columns
    return (labels_0indexed + 1).astype(np.int32)  # 1-indexed for Octave


# ── Streaming variant ──────────────────────────────────────────────────────────

_PROGRESS_RE = re.compile(r"^Harmonizing sample\s+(\d+)\s+out of\s+(\d+)")


def _is_progress_line(line: str) -> bool:
    return _PROGRESS_RE.match(line) is not None


def _parse_progress_line(line: str) -> tuple[int, int]:
    m = _PROGRESS_RE.match(line)
    if m is None:
        raise ValueError(f"Not a progress line: {line!r}")
    return int(m.group(1)), int(m.group(2))


def _write_stdin(proc: subprocess.Popen[str], data: str) -> None:
    if proc.stdin is None:
        return
    try:
        proc.stdin.write(data)
        proc.stdin.close()
    except BrokenPipeError:
        pass


def run_octave_normalize_streaming(
    pool_df: pd.DataFrame,
    nh: int,
    np_: int,
    k: int,
    octave_scripts_dir: str,
    octave_bin: str = "octave",
    timeout_s: int = 6000,
    random_seed: int | None = None,
    progress_queue: _QueueLike | None = None,
    worker_idx: int = 0,
    fixed_clusters: np.ndarray | None = None,
    skip_qn: bool = False,
) -> pd.DataFrame:
    """
    Run the full Octave normalization pipeline using a streaming subprocess.

    Identical in behaviour to run_octave_normalize() but uses subprocess.Popen
    so that stdout is read line-by-line as Octave runs. Progress lines emitted
    by Shambhala2_piped.m ("Harmonizing sample X out of Y") are intercepted and
    forwarded to progress_queue as ProgressEvent objects; all other lines are
    accumulated for parsing. A stdin-writer thread prevents deadlocks when the
    input TSV is large.

    Parameters
    ----------
    pool_df : pd.DataFrame
        Merged input + P reference matrix. Shape: (n_genes, nh + np_).
    nh : int
        Number of input samples in this batch.
    np_ : int
        Number of P reference samples.
    k : int
        Number of k-means probe clusters for CuBlock.
    octave_scripts_dir : str
        Absolute path to the directory containing the .m files.
    octave_bin : str
        Path to the Octave executable. Default: 'octave'.
    timeout_s : int
        Subprocess timeout in seconds. Default: 6000.
    random_seed : int, optional
        Seeds Octave's random state for reproducible k-means.
    progress_queue : _QueueLike, optional
        Queue receiving ProgressEvent objects. None disables reporting.
    worker_idx : int
        Worker index used to label ProgressEvent objects. Default: 0.

    Returns
    -------
    pd.DataFrame
        Quantile-normalized and CuBlock-normalized expression values.
        Shape: (n_genes, nh). Row index = gene symbols.

    Raises
    ------
    OctaveExecutionError
        If Octave exits with a non-zero return code.
    OctaveOutputError
        If stdout cannot be parsed or has unexpected shape.
    """
    input_sample_names = list(pool_df.columns[:nh])
    gene_symbols = list(pool_df.index)

    tsv_string = _format_pool_as_tsv(pool_df)

    seed_preamble = (
        f"rand('state', {random_seed}); randn('state', {random_seed}); "
        if random_seed is not None
        else ""
    )

    # Approach E: write fixed_clusters to a temp file for Octave to load
    _clusters_tmpfile = None
    clusters_preamble = ""
    if skip_qn:
        script_name = "Shambhala2_piped_preqn.m"
    else:
        script_name = "Shambhala2_piped.m"
    if fixed_clusters is not None:
        _clusters_tmpfile = tempfile.NamedTemporaryFile(
            mode="w", suffix=".txt", delete=False
        )
        for label in fixed_clusters:
            _clusters_tmpfile.write(f"{int(label)}\n")
        _clusters_tmpfile.flush()
        _clusters_tmpfile.close()
        clusters_preamble = (
            f"FIXED_CLUSTERS = int32(load('{_clusters_tmpfile.name}')); "
        )
        script_name = "Shambhala2_piped_fixed.m"

    eval_preamble = (
        f"addpath('{octave_scripts_dir}'); "
        f"{seed_preamble}"
        f"NH={nh}; NP={np_}; k={k}; "
        f"{clusters_preamble}"
        f"source('{octave_scripts_dir}/{script_name}');"
    )

    try:
        proc = subprocess.Popen(
            [octave_bin, "--no-gui", "--eval", eval_preamble],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except FileNotFoundError:
        raise OctaveExecutionError(
            f"Octave executable not found: {octave_bin!r}. "
            "Install Octave or pass the correct path via --octave-bin."
        )

    stdin_thread = threading.Thread(target=_write_stdin, args=(proc, tsv_string), daemon=True)
    stdin_thread.start()

    stdout_lines: list[str] = []
    t_start = time.time()

    for raw_line in proc.stdout:  # type: ignore[union-attr]
        line = raw_line.rstrip("\n")
        if _is_progress_line(line):
            x, y = _parse_progress_line(line)
            if progress_queue is not None:
                progress_queue.put(
                    ProgressEvent(
                        worker_idx=worker_idx,
                        sample_done=x,
                        sample_total=y,
                        elapsed_s=time.time() - t_start,
                        is_final=False,
                    )
                )
        else:
            stdout_lines.append(line)

    stdin_thread.join()

    stderr_output = proc.stderr.read() if proc.stderr is not None else ""  # type: ignore[union-attr]
    try:
        proc.wait(timeout=timeout_s)
    except subprocess.TimeoutExpired:
        proc.kill()
        raise OctaveExecutionError(
            f"Octave subprocess timed out after {timeout_s}s."
        )

    if _clusters_tmpfile is not None:
        import os as _os  # noqa: PLC0415
        try:
            _os.unlink(_clusters_tmpfile.name)
        except OSError:
            pass

    if proc.returncode != 0:
        raise OctaveExecutionError(
            f"Octave exited with code {proc.returncode}.\n"
            f"--- STDERR (full) ---\n{stderr_output}\n"
            f"--- STDOUT (first 2000 chars) ---\n{chr(10).join(stdout_lines)[:2000]}"
        )

    result = _parse_octave_stdout("\n".join(stdout_lines), gene_symbols, input_sample_names)

    if progress_queue is not None:
        progress_queue.put(
            ProgressEvent(
                worker_idx=worker_idx,
                sample_done=nh,
                sample_total=nh,
                elapsed_s=time.time() - t_start,
                is_final=True,
            )
        )

    return result
