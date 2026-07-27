"""
Unit tests for shambhala/io_utils.py.
"""

import io
import gzip
import pathlib
import tempfile
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest

import sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from shambhala.io_utils import read_expression, write_expression, _parse_s3_uri

FIXTURES_DIR = pathlib.Path(__file__).parent / "fixtures"


def _make_small_df() -> pd.DataFrame:
    rng = np.random.default_rng(0)
    data = rng.random((5, 10))
    return pd.DataFrame(data, index=[f"S{i}" for i in range(5)], columns=[f"G{j}" for j in range(10)])


def test_read_local_csv():
    df = read_expression(str(FIXTURES_DIR / "input_10samples.csv"))
    assert df.shape == (10, 34592)
    assert df.index[0] == "СTO_1"
    assert df.dtypes.iloc[0].kind == "f"


def test_read_local_tsv_gz():
    small = _make_small_df()
    with tempfile.NamedTemporaryFile(suffix=".tsv.gz", delete=False) as tf:
        path = tf.name
    try:
        # Write as gzipped TSV manually
        buf = io.StringIO()
        small.to_csv(buf, sep="\t")
        with gzip.open(path, "wt", encoding="utf-8") as gf:
            gf.write(buf.getvalue())
        result = read_expression(path)
        pd.testing.assert_frame_equal(result, small, check_dtype=False)
    finally:
        pathlib.Path(path).unlink(missing_ok=True)


def test_write_and_read_roundtrip():
    original = _make_small_df()
    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as tf:
        path = tf.name
    try:
        write_expression(original, path)
        result = read_expression(path)
        pd.testing.assert_frame_equal(result, original, check_dtype=False)
    finally:
        pathlib.Path(path).unlink(missing_ok=True)


def test_parse_s3_uri():
    bucket, key = _parse_s3_uri("s3://my-bucket/path/to/file.csv")
    assert bucket == "my-bucket"
    assert key == "path/to/file.csv"


def test_s3_read_mocked():
    small = _make_small_df()
    csv_bytes = small.to_csv().encode("utf-8")

    mock_client = MagicMock()
    mock_client.get_object.return_value = {"Body": io.BytesIO(csv_bytes)}

    with patch("shambhala.io_utils.get_s3_client", return_value=mock_client):
        result = read_expression("s3://test-bucket/data/matrix.csv")

    mock_client.get_object.assert_called_once_with(Bucket="test-bucket", Key="data/matrix.csv")
    pd.testing.assert_frame_equal(result, small, check_dtype=False)
