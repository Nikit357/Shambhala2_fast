"""
Shared constants, S3 helpers, and metric utilities for the Shambhala cross-product benchmark.

This module must NOT import rpy2. All functions are pure Python + pandas/numpy.
Metric and post-removal functions are copied verbatim from bench_shared.py /
run_norm_job.py in harmonization-scripts/ to avoid that module's rpy2 top-level import.
"""
from __future__ import annotations

import gzip
import io
import os

import botocore.exceptions
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# ── Constants ──────────────────────────────────────────────────────────────────
# Deployment-specific. Set before running any benchmark script:
#   export SHAMBHALA_S3_BUCKET=your-bucket
# No default is provided on purpose: a wrong-bucket default would fail silently
# or write to somewhere unintended, whereas a missing variable fails at import.
S3_BUCKET = os.environ["SHAMBHALA_S3_BUCKET"]
S3_PREFIX = "FL_batch_correction"
S3_CALIB_PREFIX = "FL_batch_correction/calibration_datasets"
BATCH_COL = "RNA_BATCH"
BIO_COL = "Diagnosis_cell_type_unified"

# ── Calibration dataset registry ───────────────────────────────────────────────
P_DATASETS: dict[str, str] = {
    "P0std": "P0_standard.csv",
    "ANTE": "ANTE.csv",
    "GTExAffy": "GTExAffymetrix.csv",
    "NBBags": "Normal_B_BG_BAGS.csv",
    "NBGPL570": "Normal_B_GPL570.csv",
    "NBKass": "Normal_B_Kassandra.csv",
    "NBRNAseq": "Normal_B_RNASeq.csv",
    "NBext": "Normal_B_cells_extended.csv",
    "Oncobox": "OncoboxCancer.csv",
}

Q_DATASETS: dict[str, str] = {
    "Q0std": "Q0_standard.csv",
    "QNBKass": "Normal_B_Kassandra.csv",
}

# 18 entries: all (p_key, q_key) combinations
SHAMBHALA_VARIANTS: dict[str, tuple[str, str]] = {
    f"shambhala_{p}_{q}": (p, q)
    for q in Q_DATASETS
    for p in P_DATASETS
}

ALL_STRATEGIES: list[str] = [
    "S0_no_removal",
    "A_confirmed_bad",
    "B_extended_bad",
    "C_rnaseq_only",
    "D_malignant_only",
    "E1_iterative_r1",
    "E2_iterative_r2",
    "E3_iterative_r3",
    "F_microarray_only",
    "G_affymetrix_only",
    "H_affymetrix_extended",
    "I_rare_batches_removed",
    "J_ff_only",
    "K_ffpe_only",
]
ALL_IMPUTATION: list[str] = ["strict", "knn", "softimpute"]
ALL_METHODS: list[str] = sorted(SHAMBHALA_VARIANTS.keys())


# ── S3 key functions ───────────────────────────────────────────────────────────

def s3_key_prepared_exp(strat: str, imp: str) -> str:
    return f"{S3_PREFIX}/prepared/{strat}__{imp}__exp.tsv.gz"


def s3_key_prepared_ann(strat: str, imp: str) -> str:
    return f"{S3_PREFIX}/prepared/{strat}__{imp}__ann.tsv.gz"


def s3_key_exp(strat: str, imp: str, method: str, post_rm: bool) -> str:
    suffix = "post1" if post_rm else "post0"
    return f"{S3_PREFIX}/exp/{strat}__{imp}__{method}__{suffix}.tsv.gz"


def s3_key_sidecar(imp: str, method: str) -> str:
    return f"{S3_PREFIX}/shambhala_metrics/{imp}__{method}.json"


def s3_key_calib(filename: str) -> str:
    return f"{S3_CALIB_PREFIX}/{filename}"


def s3_key_shambhala_metrics_csv() -> str:
    return f"{S3_PREFIX}/shambhala_metrics.csv"


# ── S3 helpers ─────────────────────────────────────────────────────────────────

def s3_exists(s3_client: object, key: str) -> bool:
    """Check whether an S3 key exists without downloading its content."""
    try:
        s3_client.head_object(Bucket=S3_BUCKET, Key=key)  # type: ignore[attr-defined]
        return True
    except botocore.exceptions.ClientError:
        return False


def download_exp_from_s3(s3_client: object, key: str) -> pd.DataFrame:
    """Download a gzip-compressed TSV expression matrix from S3 (samples × genes)."""
    buf = io.BytesIO()
    s3_client.download_fileobj(S3_BUCKET, key, buf)  # type: ignore[attr-defined]
    buf.seek(0)
    with gzip.open(buf, "rb") as gz:
        return pd.read_csv(gz, sep="\t", index_col=0)


def download_ann_from_s3(s3_client: object, key: str) -> pd.DataFrame:
    """Download a gzip-compressed TSV annotation DataFrame from S3."""
    buf = io.BytesIO()
    s3_client.download_fileobj(S3_BUCKET, key, buf)  # type: ignore[attr-defined]
    buf.seek(0)
    with gzip.open(buf, "rb") as gz:
        return pd.read_csv(gz, sep="\t", index_col=0, low_memory=False)


def download_calib_df(s3_client: object, filename: str) -> pd.DataFrame:
    """
    Download a calibration dataset CSV from S3.

    All calibration files use FL convention: rows = samples, columns = HGNC gene symbols.
    No transposition is needed.

    Parameters
    ----------
    s3_client
        Boto3 S3 client.
    filename
        Filename relative to the calibration_datasets/ prefix (e.g. 'P0_standard.csv').

    Returns
    -------
    pd.DataFrame
        Shape (n_samples, n_genes). Index = sample IDs, columns = HGNC gene symbols.
    """
    key = s3_key_calib(filename)
    buf = io.BytesIO()
    s3_client.download_fileobj(S3_BUCKET, key, buf)  # type: ignore[attr-defined]
    buf.seek(0)
    return pd.read_csv(buf, index_col=0)


def upload_exp_to_s3(exp_df: pd.DataFrame, s3_client: object, key: str) -> None:
    """Upload a expression DataFrame as gzip-compressed TSV to S3."""
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb") as gz:
        exp_df.to_csv(gz, sep="\t")
    buf.seek(0)
    s3_client.upload_fileobj(buf, S3_BUCKET, key)  # type: ignore[attr-defined]


# ── Memory helper ──────────────────────────────────────────────────────────────

def free_gb() -> float:
    """Return available RAM in GB, reading container cgroup limit when inside K8s."""
    def _read_int(path: str) -> int | None:
        try:
            val = open(path).read().strip()
            return None if val == "max" else int(val)
        except (OSError, ValueError):
            return None

    limit = _read_int("/sys/fs/cgroup/memory.max") or _read_int(
        "/sys/fs/cgroup/memory/memory.limit_in_bytes"
    )
    usage = _read_int("/sys/fs/cgroup/memory.current") or _read_int(
        "/sys/fs/cgroup/memory/memory.usage_in_bytes"
    )
    if limit is not None and usage is not None:
        return (limit - usage) / 1e9
    try:
        import psutil
        return psutil.virtual_memory().available / 1e9
    except ImportError:
        return 999.0


# ── Batch-effect metric ────────────────────────────────────────────────────────

def pca_variance_explained_by_batch(
    exp_df: pd.DataFrame,
    batch_series: pd.Series,
    n_components: int = 10,
) -> float:
    """
    Mean fraction of PCA variance explained by batch (one-way ANOVA R²).

    Parameters
    ----------
    exp_df
        Expression matrix (samples × genes).
    batch_series
        Batch labels indexed by sample ID.
    n_components
        Number of PCs to average over.

    Returns
    -------
    Mean R² across first n_components PCs.
    """
    pca = PCA(n_components=n_components)
    coords = pca.fit_transform(StandardScaler().fit_transform(exp_df.fillna(0)))
    groups = batch_series.reindex(exp_df.index).fillna("NA")
    r2_list = []
    for pc in range(n_components):
        vals = coords[:, pc]
        grand_mean = vals.mean()
        ss_between = sum(
            len(vals[groups == g]) * (vals[groups == g].mean() - grand_mean) ** 2
            for g in groups.unique()
            if (groups == g).sum() > 1
        )
        ss_total = ((vals - grand_mean) ** 2).sum()
        r2_list.append(ss_between / ss_total if ss_total > 0 else 0)
    return float(np.mean(r2_list))


def r2_batch(
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    batch_col: str = BATCH_COL,
    n_pcs: int = 10,
) -> float:
    """
    PCA R² for a batch column in ann_df.

    Parameters
    ----------
    exp_df
        Expression matrix (samples × genes).
    ann_df
        Annotation DataFrame containing batch_col.
    batch_col
        Column name for grouping (RNA_BATCH or Diagnosis_cell_type_unified).
    n_pcs
        Number of PCs to average over.

    Returns
    -------
    Mean R² across PCs.
    """
    return pca_variance_explained_by_batch(exp_df, ann_df[batch_col], n_pcs)


# ── Post-removal ───────────────────────────────────────────────────────────────

def identify_outlier_batches(
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    batch_col: str = BATCH_COL,
    n_pcs: int = 2,
    n_outliers: int = 1,
    min_batch_size: int = 20,
) -> list[str]:
    """
    Identify batches whose PCA centroid is furthest from the global centroid.

    Parameters
    ----------
    exp_df
        Expression matrix (samples × genes).
    ann_df
        Annotation with batch_col.
    batch_col
        Column for batch grouping.
    n_pcs
        Number of PCs for centroid calculation.
    n_outliers
        Number of outlier batches to return.
    min_batch_size
        Batches smaller than this are excluded from outlier detection.

    Returns
    -------
    List of outlier batch names, length <= n_outliers.
    """
    common = exp_df.index.intersection(ann_df.index)
    X = StandardScaler().fit_transform(exp_df.loc[common].fillna(0))
    coords = PCA(n_components=n_pcs).fit_transform(X)
    coords_df = pd.DataFrame(coords, index=common)
    batches = ann_df.loc[common, batch_col].fillna("NA")
    global_centroid = coords_df.mean(axis=0).values

    batch_centroids: dict[str, float] = {}
    for b in batches.unique():
        idx = batches[batches == b].index
        if len(idx) < min_batch_size:
            continue
        batch_centroids[b] = float(
            np.linalg.norm(coords_df.loc[idx].mean(axis=0).values - global_centroid)
        )

    ranked = sorted(batch_centroids, key=batch_centroids.get, reverse=True)  # type: ignore[arg-type]
    return ranked[:n_outliers]
