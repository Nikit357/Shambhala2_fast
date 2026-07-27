"""
Parallel sample batch dispatcher for the Shambhala2 harmonization pipeline.

Uses ProcessPoolExecutor so each worker is a fully isolated subprocess,
avoiding file-descriptor conflicts that would occur with ThreadPoolExecutor
when multiple Octave processes share stdin/stdout pipes.
"""

import logging
import multiprocessing
import pathlib
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd

from shambhala.octave_bridge import run_octave_normalize, run_octave_normalize_streaming
from shambhala.progress_display import ProgressEvent, _ProgressRenderer, _QueueLike
from shambhala import q_rescale as q_rescale_module

logger = logging.getLogger(__name__)


def _harmonize_one_batch(
    batch_sample_names: list[str],
    input_df: pd.DataFrame,
    p_df: pd.DataFrame,
    k: int,
    octave_scripts_dir: str,
    octave_bin: str,
    timeout_s: int,
    random_seed: int | None,
    fixed_clusters: np.ndarray | None = None,
    python_cublock: bool = False,
    skip_qn: bool = False,
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
    batch_df = input_df[batch_sample_names]
    pool_df = pd.concat([batch_df, p_df], axis=1)

    nh = len(batch_sample_names)
    np_ = p_df.shape[1]
    t0 = time.time()

    if python_cublock:
        result = _run_python_cublock_batch(
            pool_df=pool_df,
            batch_sample_names=batch_sample_names,
            nh=nh,
            k=k,
            random_seed=random_seed,
        )
    else:
        result = run_octave_normalize(
            pool_df=pool_df,
            nh=nh,
            np_=np_,
            k=k,
            octave_scripts_dir=octave_scripts_dir,
            octave_bin=octave_bin,
            timeout_s=timeout_s,
            random_seed=random_seed,
            fixed_clusters=fixed_clusters,
            skip_qn=skip_qn,
        )

    elapsed = time.time() - t0
    logger.info(
        "Batch done (python_cublock=%s): NH=%d, NP=%d, n_genes=%d, elapsed=%.0fs.",
        python_cublock, nh, np_, pool_df.shape[0], elapsed,
    )
    return result


def _harmonize_one_batch_with_progress(
    batch_sample_names: list[str],
    input_df: pd.DataFrame,
    p_df: pd.DataFrame,
    k: int,
    octave_scripts_dir: str,
    octave_bin: str,
    timeout_s: int,
    random_seed: int | None,
    progress_queue: _QueueLike | None,
    worker_idx: int,
    fixed_clusters: np.ndarray | None = None,
    python_cublock: bool = False,
    skip_qn: bool = False,
) -> pd.DataFrame:
    """
    Harmonize one batch of samples with per-sample progress reporting.

    Identical to _harmonize_one_batch but calls run_octave_normalize_streaming,
    which intercepts Octave's per-sample progress lines and forwards them to
    progress_queue as ProgressEvent objects.

    Must be a module-level function to be picklable by ProcessPoolExecutor.

    Parameters
    ----------
    batch_sample_names : list[str]
        Sample IDs assigned to this batch.
    input_df : pd.DataFrame
        Full input matrix (genes × samples).
    p_df : pd.DataFrame
        P reference matrix (genes × p_samples).
    k : int
        k-means clusters for CuBlock.
    octave_scripts_dir : str
        Absolute path to the octave/ directory.
    octave_bin : str
        Path to the Octave executable.
    timeout_s : int
        Per-batch Octave subprocess timeout in seconds.
    random_seed : int or None
        Passed to run_octave_normalize_streaming for reproducible k-means.
    progress_queue : _QueueLike or None
        Queue receiving ProgressEvent objects. None disables reporting.
    worker_idx : int
        Worker index used to label ProgressEvent objects.

    Returns
    -------
    pd.DataFrame
        Quantile-normalized + CuBlock-normalized batch (genes × batch_samples).
    """
    batch_df = input_df[batch_sample_names]
    pool_df = pd.concat([batch_df, p_df], axis=1)

    nh = len(batch_sample_names)
    np_ = p_df.shape[1]
    t0 = time.time()

    if python_cublock:
        result = _run_python_cublock_batch(
            pool_df=pool_df,
            batch_sample_names=batch_sample_names,
            nh=nh,
            k=k,
            random_seed=random_seed,
            progress_queue=progress_queue,
            worker_idx=worker_idx,
        )
    else:
        result = run_octave_normalize_streaming(
            pool_df=pool_df,
            nh=nh,
            np_=np_,
            k=k,
            octave_scripts_dir=octave_scripts_dir,
            octave_bin=octave_bin,
            timeout_s=timeout_s,
            random_seed=random_seed,
            progress_queue=progress_queue,
            worker_idx=worker_idx,
            fixed_clusters=fixed_clusters,
            skip_qn=skip_qn,
        )

    elapsed = time.time() - t0
    logger.info(
        "Batch done (python_cublock=%s): NH=%d, NP=%d, n_genes=%d, elapsed=%.0fs.",
        python_cublock, nh, np_, pool_df.shape[0], elapsed,
    )
    return result


def _run_python_cublock_batch(
    pool_df: pd.DataFrame,
    batch_sample_names: list[str],
    nh: int,
    k: int,
    random_seed: int | None,
    progress_queue: _QueueLike | None = None,
    worker_idx: int = 0,
) -> pd.DataFrame:
    """
    Process one batch with Python QN + Python CuBlock (Approach G).

    For each input sample, builds [sample | P] pool, applies Python QN via
    qnorm, then applies cublock_python. Emits ProgressEvent per sample if
    progress_queue is provided.

    Returns a DataFrame (genes × batch_samples) matching the Octave bridge output.
    """
    import qnorm as qnorm_lib
    from shambhala.cublock_python import cublock_python

    pool_np = pool_df.values.astype(float)
    gene_index = pool_df.index
    n_p = pool_np.shape[1] - nh
    result_cols: list[np.ndarray] = []
    t_sample_start = time.time()

    for col_idx in range(nh):
        sample_pool = np.column_stack([pool_np[:, col_idx], pool_np[:, nh : nh + n_p]])
        qn_pool = qnorm_lib.quantile_normalize(sample_pool, axis=1)
        cublock_result = cublock_python(qn_pool, n_reps=30, k=k, random_seed=random_seed)
        result_cols.append(cublock_result[:, 0])

        if progress_queue is not None:
            elapsed = time.time() - t_sample_start
            speed = elapsed / (col_idx + 1)
            progress_queue.put(
                ProgressEvent(
                    worker_idx=worker_idx,
                    done=col_idx + 1,
                    total=nh,
                    elapsed_s=elapsed,
                    speed_s_per_sample=speed,
                )
            )

    result_np = np.column_stack(result_cols) if len(result_cols) > 1 else result_cols[0][:, np.newaxis]
    return pd.DataFrame(result_np, index=gene_index, columns=batch_sample_names)


def harmonize_parallel(
    input_df: pd.DataFrame,
    p_df: pd.DataFrame,
    rm: pd.Series,
    rs: pd.Series,
    k: int = 5,
    n_workers: int = 5,
    octave_scripts_dir: str | None = None,
    octave_bin: str = "octave",
    timeout_s: int = 6000,
    random_seed: int | None = None,
    disable_progress: bool = False,
    fixed_clusters: np.ndarray | None = None,
    python_cublock: bool = False,
    skip_qn: bool = False,
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
        Must be pre-filtered to exactly the same gene set as input_df
        (i.e. p_df.index == input_df.index). Caller is responsible for
        this alignment after NA dropping.
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
    disable_progress : bool
        If True, skip Manager and renderer creation entirely. Default: False.
        Useful for unit tests where Manager startup overhead is undesirable.

    Returns
    -------
    pd.DataFrame
        Fully harmonized expression matrix. Shape: (n_genes, n_samples).
        Rows = genes, columns = samples (same order as input).

    Raises
    ------
    RuntimeError
        If any worker process raises an exception (original exception chained).
    """

    if octave_scripts_dir is None:
        octave_scripts_dir = str(
            pathlib.Path(__file__).resolve().parent.parent / "octave"
        )

    n_samples = input_df.shape[1]
    effective_workers = min(n_workers, n_samples)
    if effective_workers < n_workers:
        logger.info(
            "n_workers=%d clamped to %d (fewer samples than workers).",
            n_workers,
            effective_workers,
        )

    sample_names = list(input_df.columns)
    batches = [
        list(arr)
        for arr in np.array_split(sample_names, effective_workers)
        if len(arr) > 0
    ]

    logger.info(
        "Dispatching %d samples across %d workers (%d–%d samples/worker).",
        n_samples,
        effective_workers,
        min(len(b) for b in batches),
        max(len(b) for b in batches),
    )

    # ── Progress display setup ────────────────────────────────────────────────
    manager: multiprocessing.managers.SyncManager | None = None
    progress_queue: _QueueLike | None = None
    renderer: _ProgressRenderer | None = None

    if not disable_progress:
        manager = multiprocessing.Manager()
        raw_queue = manager.Queue()  # type: ignore[var-annotated]
        progress_queue = raw_queue
        renderer = _ProgressRenderer(
            queue=progress_queue,
            n_workers=effective_workers,
            total_samples=n_samples,
            worker_batch_sizes=[len(b) for b in batches],
            use_ansi=sys.stderr.isatty(),
        )
        renderer.daemon = True
        renderer.start()

    # ── Dispatch workers ──────────────────────────────────────────────────────
    batch_results: dict[int, pd.DataFrame] = {}

    # Record child PIDs before spawning so we can terminate only our workers on
    # unexpected exits (timeout, KeyboardInterrupt, etc.).
    _children_before = {p.pid for p in multiprocessing.active_children()}
    executor = ProcessPoolExecutor(max_workers=effective_workers)

    try:
        future_to_idx = {
            executor.submit(
                _harmonize_one_batch_with_progress,
                batch,
                input_df,
                p_df,
                k,
                octave_scripts_dir,
                octave_bin,
                timeout_s,
                random_seed,
                progress_queue,
                idx,
                fixed_clusters,
                python_cublock,
                skip_qn,
            ): idx
            for idx, batch in enumerate(batches)
        }

        for future in as_completed(future_to_idx):
            idx = future_to_idx[future]
            try:
                batch_results[idx] = future.result()
                logger.info("Worker %d/%d complete.", idx + 1, effective_workers)
            except Exception as exc:
                raise RuntimeError(
                    f"Worker for batch {idx} (samples: {batches[idx]}) failed."
                ) from exc

    except BaseException:
        # Terminate our worker processes immediately so shutdown() does not
        # block waiting for them (e.g. after pytest-timeout raises Timeout).
        for p in multiprocessing.active_children():
            if p.pid not in _children_before:
                p.terminate()
        executor.shutdown(wait=False, cancel_futures=True)
        raise

    else:
        executor.shutdown(wait=True)

    finally:
        if renderer is not None:
            renderer.stop()
        if manager is not None:
            manager.shutdown()

    # Reassemble in original batch order, then restore original sample column order
    ordered_results = [batch_results[i] for i in range(len(batches))]
    assembled = pd.concat(ordered_results, axis=1)[sample_names]

    harmonized = q_rescale_module.rescale(assembled, rm, rs)
    return harmonized
