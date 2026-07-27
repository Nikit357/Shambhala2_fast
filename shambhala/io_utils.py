"""
Unified file I/O for expression matrices.

Abstracts local filesystem and S3 object storage behind a single interface.
All read/write operations use the FL data convention:
    rows = samples (index = sample IDs)
    columns = genes (headers = HGNC gene symbols)
"""

import io
import pathlib

import pandas as pd


def _parse_s3_uri(uri: str) -> tuple[str, str]:
    """
    Parse an S3 URI into (bucket_name, key).

    Parameters
    ----------
    uri : str
        S3 URI of the form 's3://bucket-name/path/to/file.csv'.

    Returns
    -------
    tuple[str, str]
        (bucket_name, key) e.g. ('my-bucket', 'path/to/file.csv').

    Raises
    ------
    ValueError
        If the URI does not start with 's3://'.
    """
    if not uri.startswith("s3://"):
        raise ValueError(f"Not an S3 URI: {uri!r}")
    without_prefix = uri[len("s3://"):]
    bucket, _, key = without_prefix.partition("/")
    return bucket, key


def get_s3_client():
    """
    Return a boto3 S3 client using the default credential chain.

    Credentials are resolved in order: environment variables, ~/.aws/credentials,
    IAM instance role. No credentials are hardcoded.

    Returns
    -------
    boto3.client
        Configured S3 client.
    """
    import boto3
    return boto3.client("s3")


def _infer_separator(path: str) -> str:
    """
    Infer the column separator from the file extension.

    Parameters
    ----------
    path : str
        File path or S3 URI.

    Returns
    -------
    str
        ',' for .csv or .csv.gz; '\t' for .tsv or .tsv.gz.
    """
    lower = path.lower()
    if lower.endswith(".tsv") or lower.endswith(".tsv.gz"):
        return "\t"
    return ","


def read_expression(path: str, sep: str = None) -> pd.DataFrame:
    """
    Read an expression matrix from a local path or S3 URI.

    The file must follow the FL project convention:
        - Rows = samples (row index = sample IDs)
        - Columns = genes (column headers = HGNC gene symbols)
        - First column of the CSV is the sample ID (used as the index)

    Accepts .csv, .csv.gz, .tsv, .tsv.gz. Separator is inferred from extension
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
    resolved_sep = sep if sep is not None else _infer_separator(path)
    compression = "gzip" if path.lower().endswith(".gz") else None

    if path.startswith("s3://"):
        bucket, key = _parse_s3_uri(path)
        client = get_s3_client()
        response = client.get_object(Bucket=bucket, Key=key)
        raw_bytes = response["Body"].read()
        buffer = io.BytesIO(raw_bytes)
        return pd.read_csv(
            buffer, sep=resolved_sep, index_col=0,
            compression=compression, low_memory=False,
        )

    return pd.read_csv(
        path, sep=resolved_sep, index_col=0,
        compression=compression, low_memory=False,
    )


def write_expression(df: pd.DataFrame, path: str, sep: str = ",") -> None:
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
    if path.startswith("s3://"):
        bucket, key = _parse_s3_uri(path)
        client = get_s3_client()
        buffer = io.StringIO()
        df.to_csv(buffer, sep=sep)
        client.put_object(
            Bucket=bucket,
            Key=key,
            Body=buffer.getvalue().encode("utf-8"),
        )
        return

    parent = pathlib.Path(path).parent
    if not parent.exists():
        parent.mkdir(parents=True, exist_ok=True)

    df.to_csv(path, sep=sep)
