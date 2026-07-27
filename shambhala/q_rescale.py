"""
Q-rescaling step for the Shambhala2 harmonization pipeline.

Translates the final rescaling step from Shambhala2.R into pandas/numpy,
eliminating the R dependency. The algorithm is a faithful translation of
the original R rowMeans / rowSds logic.
"""

import logging
import warnings

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


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
    q_values = q_df.values.astype(float)

    if (q_values < 0).any():
        raise ValueError(
            "Q reference matrix contains negative values. "
            "Expression values must be non-negative."
        )

    gene_symbols = list(q_df.columns)

    if q_pseudocount > 0:
        log_q = np.log(q_values + q_pseudocount)
        rm = pd.Series(np.mean(log_q, axis=0), index=gene_symbols)
        rs = pd.Series(np.std(log_q, axis=0, ddof=1), index=gene_symbols)
        return rm, rs

    # q_pseudocount == 0: exclude genes with any zero or negative value
    has_zero = (q_values == 0).any(axis=0)
    n_zero = int(has_zero.sum())
    if n_zero > 0:
        excluded_genes = [gene_symbols[i] for i in range(len(gene_symbols)) if has_zero[i]]
        warnings.warn(
            f"{n_zero} gene(s) in Q have zero counts and are excluded from Q statistics "
            f"(q_pseudocount=0). Excluded: {excluded_genes[:5]}"
            + (" ..." if n_zero > 5 else ""),
            RuntimeWarning,
            stacklevel=2,
        )
        valid_mask = ~has_zero
        q_valid = q_values[:, valid_mask]
        gene_symbols_valid = [g for g, keep in zip(gene_symbols, valid_mask) if keep]
    else:
        q_valid = q_values
        gene_symbols_valid = gene_symbols

    log_q = np.log(q_valid)
    rm = pd.Series(np.mean(log_q, axis=0), index=gene_symbols_valid)
    rs = pd.Series(np.std(log_q, axis=0, ddof=1), index=gene_symbols_valid)
    return rm, rs


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
    common_genes = octave_output.index.intersection(rm.index)
    n_dropped = len(octave_output.index) - len(common_genes)
    if n_dropped > 0:
        logger.warning(
            "%d gene(s) in Octave output are absent from Q statistics and will be dropped.",
            n_dropped,
        )

    aligned_output = octave_output.loc[common_genes]
    aligned_rm = rm.loc[common_genes].values.reshape(-1, 1)
    aligned_rs = rs.loc[common_genes].values.reshape(-1, 1)

    log_octave = np.log(aligned_output.values + 1)
    rescaled_log = aligned_rm + aligned_rs * log_octave
    rescaled = np.exp(rescaled_log)

    return pd.DataFrame(rescaled, index=common_genes, columns=octave_output.columns)
