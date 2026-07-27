# Pod Python Package Error — Root Cause Analysis and Fix Plan

**Date:** 2026-05-15  
**Symptom:** `pip3 install -r requirements.txt` fails inside `shambhala-bench` with
`ModuleNotFoundError: No module named 'pkg_resources'` while building pandas.

---

## 1. What the Logs Show (Attempt 2)

Step 1 (SSH) and step 2 (Octave + Python) now succeed — the `liboctave-dev` fix from attempt 1
worked. The failure is in step 3 (pip install):

```
Collecting pandas<2.0,>=1.5.3 (from -r /etc/shambhala-deps/requirements.txt (line 3))
  Downloading pandas-1.5.3.tar.gz (5.2 MB)            ← source tarball, not a wheel
  Installing build dependencies: finished with status 'done'
  Getting requirements to build wheel: finished with status 'error'

  × Getting requirements to build wheel did not run successfully.
  │ exit code: 1
  ╰─> [20 lines of output]
      Traceback (most recent call last):
        File ".../setuptools/build_meta.py", line 317, in run_setup
          exec(code, locals())
        File "<string>", line 19, in <module>
      ModuleNotFoundError: No module named 'pkg_resources'
```

---

## 2. Root Cause #1 — pandas 1.5.3 Has No Wheel for Python 3.12

This is the primary cause. Everything else is a consequence of it.

**Timeline:**
- pandas 1.5.3 was released **January 2023** — pre-built wheels exist for Python 3.8, 3.9, 3.10, and 3.11 only.
- Python 3.12 was released **October 2023** — after pandas 1.5.x was already abandoned. No `cp312` wheel was ever published for pandas 1.5.x.
- Ubuntu 24.04 (the pod base image) ships **Python 3.12.3**.

When pip resolves `pandas>=1.5.3,<2.0` on Python 3.12, the only version in range that exists is 1.5.3, but there is no pre-built binary wheel for `cp312`. pip falls back to downloading the source tarball (`pandas-1.5.3.tar.gz`, 5.2 MB) and building it from source.

---

## 3. Root Cause #2 — `pkg_resources` Is Not Available in the Isolated Build Environment

pip uses **build isolation** when compiling from source (PEP 517). It creates a temporary virtual
environment (`/tmp/pip-build-env-.../`) and installs the build dependencies listed in
`pyproject.toml`. pandas 1.5.3's `setup.py` then runs inside this isolated subprocess.

The `setup.py` for pandas 1.5.3 (line 19) contains:
```python
from pkg_resources import parse_version
```

`pkg_resources` is provided by `setuptools`. However, the isolated build environment's subprocess
does not have `pkg_resources` available in its Python path for one of these reasons:
- Modern pip (24.0) uses a stricter build isolation sandbox where `setuptools` is installed but
  the subprocess spawned by `setuptools/build_meta.py` cannot find `pkg_resources` in its
  modified path.
- The system pip from `/usr/lib/python3/dist-packages/pip` (Ubuntu's APT-managed pip) has
  specific sandboxing that conflicts with setuptools' subprocess mechanism.

**This error is a consequence of root cause #1.** If pandas had a pre-built wheel, no source
build would be attempted and `pkg_resources` would never be imported.

---

## 4. Why pandas 2.x Fixes Both Issues

pandas 2.x (released April 2023) added Python 3.12 wheel support in 2.1.0 and has full support
in 2.2.x. When pip finds a pre-built wheel:
1. No source compilation happens.
2. No `pkg_resources` is imported during install.
3. Install completes in seconds instead of failing.

---

## 5. Compatibility Check: Does Our Code Work With pandas 2.x?

The shambhala pipeline uses only standard pandas public APIs. Checked all `.py` files in
`shambhala/` and `harmonization_scripts/`:

| Breaking change in pandas 2.0 | Used in our code? |
|---|---|
| `DataFrame.append()` removed | **No** — all `.append()` calls are on Python lists |
| `DataFrame.iteritems()` removed | No |
| `Series.iteritems()` removed | No |
| `convert_objects()` removed | No |
| Copy-on-Write (CoW) semantics changed | No explicit in-place mutations |
| `Index` integer indexing changed | No integer index lookups |

APIs we use — all stable across 1.5 → 2.x:
- `pd.read_csv()`, `pd.concat()`, `pd.DataFrame()`, `pd.Series()`
- `.loc[]`, `.index`, `.columns`, `.shape`
- `.index.intersection()`, `.isnull()`, `.fillna()`, `.T`
- `np.mean()`, `np.std()` on DataFrame values (numpy-side, unaffected)

**No code changes required** when upgrading pandas to 2.x.

---

## 6. Secondary Concern: `awscli>=1.29` Version Conflict

`awscli 1.x` pins boto3 and botocore to specific older versions. With `boto3>=1.28.17` and
`botocore>=1.31.17` unconstrained above, pip may resolve a boto3 version that awscli 1.x cannot
accept, producing a conflict. This did not surface in the current logs (the install failed before
resolving aws packages), but it may appear in attempt 3.

`awscli` is not used by any Python code in the pipeline — it's only needed for shell-level
`aws s3 ls` monitoring commands. The shell `aws` command can be installed separately via
`apt-get install awscli` (Ubuntu 24.04 ships awscli 1.x from apt, pinned to compatible boto3).

---

## 7. Proposed Fix

### Change 1 — Update pandas version in the ConfigMap

In `harmonization_scripts/k8s/shambhala-pod.yaml`, ConfigMap `shambhala-deps`:

**Before:**
```yaml
  requirements.txt: |
    boto3>=1.28.17
    botocore>=1.31.17
    pandas>=1.5.3,<2.0
    numpy>=1.24,<2.0
    scikit-learn>=1.3
    psutil>=5.9
    awscli>=1.29
```

**After:**
```yaml
  requirements.txt: |
    boto3>=1.28.17
    botocore>=1.31.17
    pandas>=2.2,<3.0
    numpy>=1.24,<2.0
    scikit-learn>=1.3
    psutil>=5.9
```

Changes:
- `pandas>=1.5.3,<2.0` → `pandas>=2.2,<3.0` — Python 3.12 wheel exists, no source build needed.
- `awscli>=1.29` — **removed** from Python requirements. Install via apt in the startup script
  instead (`apt-get install -y -qq awscli`), which avoids the boto3/botocore version conflict.
  The `aws` CLI is only needed for shell monitoring commands, not for Python imports.

### Change 2 — Add awscli to the apt-get command in step 2

```bash
apt-get install -y -qq build-essential python3 python3-pip python3-dev python-is-python3 \
    octave octave-statistics awscli
```

`awscli` from Ubuntu 24.04's apt is pre-pinned to compatible botocore/boto3 versions and does
not pollute the Python package environment.

---

## 8. Version Compatibility Summary

| Package | Old constraint | New constraint | Reason for change |
|---|---|---|---|
| `pandas` | `>=1.5.3,<2.0` | `>=2.2,<3.0` | No cp312 wheel for 1.5.x |
| `numpy` | `>=1.24,<2.0` | unchanged | numpy 1.26.x has cp312 wheel |
| `scikit-learn` | `>=1.3` | unchanged | 1.3+ has cp312 wheel |
| `boto3` | `>=1.28.17` | unchanged | pure Python, no wheel issue |
| `botocore` | `>=1.31.17` | unchanged | pure Python, no wheel issue |
| `psutil` | `>=5.9` | unchanged | has cp312 wheel |
| `awscli` | `>=1.29` (pip) | removed from pip; add to apt | avoids boto3 version conflict |

---

## 9. Expected Outcome After Fix

Step 3 log should show:
```
Collecting pandas<3.0,>=2.2 (from -r /etc/shambhala-deps/requirements.txt (line 3))
  Downloading pandas-2.2.x-cp312-cp312-linux_x86_64.whl (...)   ← pre-built wheel
...
Successfully installed boto3-... botocore-... numpy-... pandas-2.2.x scikit-learn-... psutil-...
[startup] Python imports OK
[startup] Python packages installed
[startup] Environment ready. Sync code with rsync and run.
```

No source builds, no `pkg_resources`, install completes in under 2 minutes.
