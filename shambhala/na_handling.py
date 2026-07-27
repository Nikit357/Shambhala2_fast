"""
NA detection and handling for gene expression matrices.

The CuBlock algorithm explicitly requires NaN-free input. This module provides
two strategies to satisfy that requirement before passing data to Octave.
"""

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from sklearn.impute import KNNImputer


@dataclass
class NAMask:
    """
    Metadata produced by apply_na_strategy(); needed to reconstruct the original gene set.

    Fields
    ------
    strategy : str
        'drop' or 'knn' — the strategy that was applied.
    dropped_genes : list[str]
        Gene names removed from the matrix (present for both strategies;
        for 'knn' these are genes that exceeded max_na_frac).
    imputed_positions : pd.DataFrame
        Boolean mask of imputed cell positions (strategy='knn' only).
        Same index and columns as the cleaned DataFrame.
        Empty DataFrame for strategy='drop'.
    original_gene_order : list[str]
        All gene column names in the original column order (before any dropping).
        Used by restore_na_genes() to reconstruct the full matrix.
    """

    strategy: str
    dropped_genes: list[str]
    imputed_positions: pd.DataFrame
    original_gene_order: list[str]


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
    return list(df.columns[df.isna().any(axis=0)])


def apply_na_strategy(
    df: pd.DataFrame,
    strategy: str = "drop",
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

    Raises
    ------
    ValueError
        If strategy is not 'drop' or 'knn'.
    """
    if strategy not in ("drop", "knn"):
        raise ValueError(f"Unknown NA strategy: {strategy!r}. Must be 'drop' or 'knn'.")

    original_gene_order = list(df.columns)

    if strategy == "drop":
        na_genes = detect_na_genes(df)
        clean_df = df.drop(columns=na_genes)
        mask = NAMask(
            strategy="drop",
            dropped_genes=na_genes,
            imputed_positions=pd.DataFrame(),
            original_gene_order=original_gene_order,
        )
        return clean_df, mask

    # strategy == 'knn'
    na_frac = df.isna().mean(axis=0)
    high_na_genes = list(na_frac[na_frac > max_na_frac].index)
    df_filtered = df.drop(columns=high_na_genes)

    # Record positions that will be imputed (have remaining NaN values)
    imputed_positions = df_filtered.isna().copy()

    imputer = KNNImputer(n_neighbors=knn_k)
    imputed_values = imputer.fit_transform(df_filtered.values)
    clean_df = pd.DataFrame(
        imputed_values,
        index=df_filtered.index,
        columns=df_filtered.columns,
    )

    mask = NAMask(
        strategy="knn",
        dropped_genes=high_na_genes,
        imputed_positions=imputed_positions,
        original_gene_order=original_gene_order,
    )
    return clean_df, mask


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
    if mask.strategy == "knn":
        return harmonized_df

    if not mask.dropped_genes:
        return harmonized_df[mask.original_gene_order]

    nan_block = pd.DataFrame(
        np.nan,
        index=harmonized_df.index,
        columns=mask.dropped_genes,
    )
    full_df = pd.concat([harmonized_df, nan_block], axis=1)
    return full_df[mask.original_gene_order]
