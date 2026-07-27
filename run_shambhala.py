"""
run_shambhala.py — CLI entry point for the Shambhala2 harmonization pipeline.

Harmonizes a gene expression matrix using the Shambhala2 algorithm:
    1. Quantile normalization (each sample jointly with the P calibration reference)
    2. CuBlock normalization (k-means gene clustering + cubic polynomial fitting)
    3. Q-rescaling (per-gene mean/std from the Q reference dataset)

All I/O uses the FL data convention: rows = samples, columns = genes.
S3 URIs (s3://bucket/path) are accepted for all file arguments.

Usage:
    python run_shambhala.py --input INPUT --P P_REF --Q Q_REF --output OUTPUT [options]
Exemplar real world usage:
    python run_shambhala.py --input s3://your-bucket/FL_batch_correction/prepared/S0_no_removal__strict__exp.tsv.gz --P s3://your-bucket/FL_batch_correction/calibration_datasets/P0_standard.csv --Q s3://your-bucket/FL_batch_correction/calibration_datasets/Q0_standard.csv --output s3://your-bucket/FL_batch_correction/exp/S0_no_removal__strict__shambhala_P0std_Q0std_speed_up__post0.tsv.gz --n-workers 4 --precompute-qn-reference --precompute-cublock-clusters
"""

import argparse
import logging
import pathlib
import sys
import time

import numpy as np
import pandas as pd

# Allow running as a script from within the Shambhala_containerized/ directory
_THIS_DIR = pathlib.Path(__file__).resolve().parent
if str(_THIS_DIR) not in sys.path:
    sys.path.insert(0, str(_THIS_DIR))

from shambhala import io_utils, na_handling, q_rescale, parallel as parallel_module

OCTAVE_DIR = _THIS_DIR / "octave"


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_shambhala.py",
        description=(
            "Harmonize a gene expression matrix using the Shambhala2 algorithm "
            "(quantile normalization + CuBlock + Q-rescaling, one sample at a time)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "--input",
        required=True,
        metavar="PATH",
        help=(
            "Input expression matrix. Local path or S3 URI. "
            "Convention: rows=samples, columns=genes. "
            "Formats: .csv, .tsv, .tsv.gz (auto-detected)."
        ),
    )
    required.add_argument(
        "--P",
        required=True,
        metavar="PATH",
        help=(
            "Calibration reference dataset P. Local path or S3 URI. "
            "P is used as the quantile-normalization anchor: each input "
            "sample is jointly quantile-normalized with all P samples."
        ),
    )
    required.add_argument(
        "--Q",
        required=True,
        metavar="PATH",
        help=(
            "Definitive reference dataset Q. Local path or S3 URI. "
            "Q defines the universal expression shape: per-gene mean "
            "and std of log(Q + pseudocount) are used to rescale output."
        ),
    )
    required.add_argument(
        "--output",
        required=True,
        metavar="PATH",
        help="Destination for the harmonized output matrix. Local path or S3 URI. Written in CSV format.",
    )

    optional = parser.add_argument_group("optional arguments")
    optional.add_argument(
        "--k",
        type=int,
        default=5,
        metavar="INT",
        help=(
            "Number of k-means gene clusters for the CuBlock algorithm. "
            "Default: 5. Do not change without benchmarking."
        ),
    )
    optional.add_argument(
        "--n-workers",
        type=int,
        default=5,
        metavar="INT",
        dest="n_workers",
        help=(
            "Number of parallel worker processes. Each worker runs an "
            "independent Octave subprocess on a separate batch of samples. "
            "Default: 5. Range: [1, 10]."
        ),
    )
    optional.add_argument(
        "--na-strategy",
        choices=["drop", "knn"],
        default="drop",
        dest="na_strategy",
        help=(
            "Strategy for handling NaN values in the input matrix. "
            "drop (default): Remove genes with any NaN; restored as NaN in output. "
            "knn: Impute using KNN (sklearn KNNImputer); imputed values persist."
        ),
    )
    optional.add_argument(
        "--knn-k",
        type=int,
        default=5,
        metavar="INT",
        dest="knn_k",
        help="Number of neighbours for KNN imputation. Only used with --na-strategy knn. Default: 5.",
    )
    optional.add_argument(
        "--max-na-frac",
        type=float,
        default=0.20,
        metavar="FLOAT",
        dest="max_na_frac",
        help=(
            "Maximum allowed NaN fraction per gene before dropping. "
            "Only used with --na-strategy knn. Default: 0.20."
        ),
    )
    optional.add_argument(
        "--q-pseudocount",
        type=float,
        default=1e-6,
        metavar="FLOAT",
        dest="q_pseudocount",
        help=(
            "Small constant added to Q values before log to prevent log(0)=-Inf. "
            "Default: 1e-6. Set to 0 to disable (zero-count genes in Q are then excluded)."
        ),
    )
    optional.add_argument(
        "--octave-bin",
        default="octave",
        metavar="PATH",
        dest="octave_bin",
        help=(
            "Path to the Octave executable. Default: 'octave' (resolved via PATH). "
            "For K8s pods: /opt/conda/bin/octave."
        ),
    )
    optional.add_argument(
        "--timeout-s",
        type=int,
        default=6000,
        metavar="INT",
        dest="timeout_s",
        help="Per-batch Octave subprocess timeout in seconds. Default: 6000.",
    )
    optional.add_argument(
        "--random-seed",
        type=int,
        default=None,
        metavar="INT",
        dest="random_seed",
        help=(
            "Integer seed for Octave's random state (rand and randn). "
            "Ensures reproducible k-means initialization in CuBlock. "
            "Default: not set (random initialization)."
        ),
    )
    optional.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING"],
        default="INFO",
        dest="log_level",
        help="Logging verbosity. Default: INFO.",
    )

    speedup = parser.add_argument_group("speed-up options")
    speedup.add_argument(
        "--precompute-qn-reference",
        action="store_true",
        default=False,
        dest="precompute_qn_reference",
        help=(
            "Precompute the QN reference distribution from P once in Python (qnorm library). "
            "Each sample is QN-normalized Python-side before the Octave call. "
            "When combined with --synthetic-cublock-p, Octave sees a G×2 matrix instead of G×(1+NP). "
            "~10-20× speedup for large P (e.g. NBGPL570). Near-lossless: sample contributes 1/NP to ref."
        ),
    )
    speedup.add_argument(
        "--synthetic-cublock-p",
        action="store_true",
        default=False,
        dest="synthetic_cublock_p",
        help=(
            "Pass a single synthetic P reference column to Octave CuBlock instead of all NP P samples. "
            "Synthetic P = mean of sorted P columns per gene rank. "
            "Best used with --precompute-qn-reference for maximum speedup."
        ),
    )
    speedup.add_argument(
        "--max-p-samples",
        type=int,
        default=None,
        metavar="N",
        dest="max_p_samples",
        help=(
            "Subsample P to N centroid-closest samples before harmonization. "
            "Default: disabled (full P used). "
            "Speedup proportional to NP reduction. Moderate authenticity impact."
        ),
    )
    speedup.add_argument(
        "--precompute-cublock-clusters",
        action="store_true",
        default=False,
        dest="precompute_cublock_clusters",
        help=(
            "Precompute k-means gene cluster assignments from P once; reuse per sample (CuBlock_fixed.m). "
            "Eliminates 30 k-means runs per sample. ~5-10× CuBlock speedup. "
            "Moderate authenticity impact: no per-sample k-means averaging."
        ),
    )
    speedup.add_argument(
        "--python-cublock",
        action="store_true",
        default=False,
        dest="python_cublock",
        help=(
            "Use Python (sklearn + numpy) CuBlock instead of Octave. "
            "Eliminates Octave subprocess overhead. ~3-10× speedup. "
            "Requires validation: sklearn KMeans uses k-means++ vs random init."
        ),
    )

    return parser


def _configure_logging(level_name: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level_name),
        format="[%(levelname)s] %(message)s",
        stream=sys.stdout,
    )


def main() -> None:
    parser = _build_parser()
    args = parser.parse_args()

    _configure_logging(args.log_level)
    logger = logging.getLogger(__name__)

    # Step 1 — Validate arguments
    if not (1 <= args.n_workers <= 10):
        parser.error(f"--n-workers must be in [1, 10], got {args.n_workers}.")
    if args.python_cublock and not args.precompute_qn_reference:
        logger.warning(
            "--python-cublock without --precompute-qn-reference: QN is applied inside the "
            "Python CuBlock worker (per-sample QN with P as reference)."
        )

    wall_start = time.perf_counter()

    # Step 2 — Read input matrices
    logger.info("Step 2/11 — Reading input: %s", args.input)
    input_df = io_utils.read_expression(args.input)
    logger.info("  Input shape: %d samples × %d genes", *input_df.shape)

    logger.info("Step 2/11 — Reading P reference: %s", args.P)
    p_df = io_utils.read_expression(args.P)
    p_df = p_df.astype(float)
    logger.info("  P shape: %d samples × %d genes", *p_df.shape)

    logger.info("Step 2/11 — Reading Q reference: %s", args.Q)
    q_df = io_utils.read_expression(args.Q)
    q_df = q_df.astype(float)
    logger.info("  Q shape: %d samples × %d genes", *q_df.shape)

    if p_df.isna().any().any():
        n_na_cols = int(p_df.isna().any(axis=0).sum())
        logger.error(
            "P calibration dataset has NaN values in %d gene columns. Fix the CSV.", n_na_cols
        )
        sys.exit(1)
    if q_df.isna().any().any():
        n_na_cols = int(q_df.isna().any(axis=0).sum())
        logger.error(
            "Q calibration dataset has NaN values in %d gene columns. Fix the CSV.", n_na_cols
        )
        sys.exit(1)

    # Step 3 — Transpose to internal convention (genes × samples)
    logger.info("Step 3/11 — Transposing to internal convention (genes × samples).")
    input_T = input_df.T
    p_T = p_df.T
    q_T = q_df.T

    # Step 4 — Gene intersection across input, P, Q
    logger.info("Step 4/11 — Computing gene intersection across input, P, and Q.")
    common_genes = input_T.index.intersection(p_T.index).intersection(q_T.index)
    n_input_only = len(input_T.index) - len(input_T.index.intersection(common_genes))
    n_p_only = len(p_T.index) - len(p_T.index.intersection(common_genes))
    n_q_only = len(q_T.index) - len(q_T.index.intersection(common_genes))
    logger.info(
        "  Common genes: %d (dropped %d from input, %d from P, %d from Q).",
        len(common_genes), n_input_only, n_p_only, n_q_only,
    )
    if len(common_genes) == 0:
        logger.error(
            "Gene intersection is empty. Input, P, and Q must share gene names. "
            "Check that all CSVs use HGNC gene symbols in FL convention (rows=samples, cols=genes)."
        )
        sys.exit(1)

    input_T = input_T.loc[common_genes]
    p_T = p_T.loc[common_genes]
    q_T = q_T.loc[common_genes]

    # Step 4b — Approach C: centroid-based P subsampling
    if args.max_p_samples is not None and p_T.shape[1] > args.max_p_samples:
        original_p_size = p_T.shape[1]
        p_df_aligned = p_T.T
        p_mean = p_df_aligned.mean(axis=0)
        distances = ((p_df_aligned - p_mean) ** 2).sum(axis=1)
        selected_idx = distances.nsmallest(args.max_p_samples).index
        p_T = p_T[selected_idx]
        logger.info(
            "Step 4b/11 — P subsampled from %d to %d centroid-closest samples (--max-p-samples).",
            original_p_size,
            args.max_p_samples,
        )

    # Step 5 — Apply NA strategy to input
    logger.info(
        "Step 5/11 — Applying NA strategy '%s' to input matrix.",
        args.na_strategy,
    )
    clean_df, na_mask = na_handling.apply_na_strategy(
        input_T.T,
        strategy=args.na_strategy,
        knn_k=args.knn_k,
        max_na_frac=args.max_na_frac,
    )
    clean_T = clean_df.T
    p_T = p_T.loc[clean_T.index]
    logger.info(
        "  Genes affected: %d dropped, strategy='%s'.",
        len(na_mask.dropped_genes),
        na_mask.strategy,
    )

    # Step 5b — Approaches A + D: precompute QN reference from P; apply Python-side QN
    p_for_octave = p_T
    if args.precompute_qn_reference or args.synthetic_cublock_p:
        p_sorted = np.sort(p_T.values, axis=0)   # (G, NP) — sort each P column (genes as rows)
        p_qn_reference = p_sorted.mean(axis=1)    # (G,) — mean expression per rank position

    if args.precompute_qn_reference:
        import qnorm as _qnorm  # noqa: PLC0415

        # QN each input sample independently against the fixed P reference distribution.
        # axis=1: normalize each column (sample); target: fixed sorted reference distribution.
        clean_T = _qnorm.quantile_normalize(clean_T, axis=1, target=p_qn_reference)
        logger.info(
            "Step 5b/11 — Precomputed QN reference from P (%d samples); Python-side QN applied.",
            p_T.shape[1],
        )

    if args.synthetic_cublock_p:
        p_for_octave = pd.DataFrame(
            p_qn_reference[:, np.newaxis],
            index=p_T.index,
            columns=["P_qn_reference"],
        )
        logger.info(
            "Step 5b/11 — Synthetic 1-column P passed to Octave CuBlock (--synthetic-cublock-p)."
        )

    # Step 5c — Approach E: precompute k-means gene cluster assignments from P
    fixed_clusters = None
    if args.precompute_cublock_clusters:
        from shambhala.octave_bridge import precompute_clusters_from_p  # noqa: PLC0415

        logger.info(
            "Step 5c/11 — Precomputing k-means gene clusters from P (%d samples, k=%d) ...",
            p_for_octave.shape[1],
            args.k,
        )
        fixed_clusters = precompute_clusters_from_p(
            p_T=p_for_octave,
            k=args.k,
            random_seed=args.random_seed,
        )
        logger.info("  Gene cluster assignments precomputed (k=%d).", args.k)

    # Step 6 — Pre-compute Q statistics
    logger.info("Step 6/11 — Computing Q statistics (q_pseudocount=%.2e).", args.q_pseudocount)
    q_df_filtered = q_T.loc[clean_T.index].T
    rm, rs = q_rescale.compute_q_statistics(q_df_filtered, q_pseudocount=args.q_pseudocount)
    logger.info("  Q statistics computed for %d genes.", len(rm))

    # Step 7 — Dispatch parallel harmonization
    logger.info(
        "Step 7/11 — Dispatching %d samples across %d workers.",
        clean_T.shape[1],
        args.n_workers,
    )
    harmonized_T = parallel_module.harmonize_parallel(
        input_df=clean_T,
        p_df=p_for_octave,
        rm=rm,
        rs=rs,
        k=args.k,
        n_workers=args.n_workers,
        octave_scripts_dir=str(OCTAVE_DIR),
        octave_bin=args.octave_bin,
        timeout_s=args.timeout_s,
        random_seed=args.random_seed,
        fixed_clusters=fixed_clusters,
        python_cublock=args.python_cublock,
        skip_qn=args.precompute_qn_reference,
    )
    logger.info("  Harmonization complete. Shape: %d genes × %d samples.", *harmonized_T.shape)

    # Step 8 — Restore NA genes and transpose back to FL convention (samples × genes)
    logger.info("Step 8/11 — Restoring NA genes (strategy='%s').", na_mask.strategy)
    result_df = na_handling.restore_na_genes(harmonized_T.T, na_mask)

    # Step 10 — Write output
    logger.info("Step 10/11 — Writing output to: %s", args.output)
    io_utils.write_expression(result_df, args.output)

    # Step 11 — Summary
    wall_elapsed = time.perf_counter() - wall_start
    logger.info("Step 11/11 — Done.")
    logger.info("  Total samples processed: %d", result_df.shape[0])
    logger.info("  Total genes in output: %d", result_df.shape[1])
    logger.info("  NA strategy: %s (%d genes affected)", na_mask.strategy, len(na_mask.dropped_genes))
    logger.info("  Workers used: %d", args.n_workers)
    logger.info("  Wall-clock time: %.1f s", wall_elapsed)


if __name__ == "__main__":
    main()
