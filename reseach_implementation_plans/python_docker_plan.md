# Implementation Plan: Switch Base Image to python:3.12-bookworm

**Feature:** Replace `ubuntu:24.04` base image with the official Python 3.12 image in `harmonization_scripts/k8s/shambhala-pod.yaml`.

**Goal:** Reduce startup time and installation complexity by eliminating redundant package installs that come pre-bundled in the Python image.

---

## Selected Image

**`python:3.12-bookworm`**

### Justification

The official Python registry offers four Debian-based 3.12 variants:

| Image | Base | Pre-installed packages | Size |
|---|---|---|---|
| `python:3.12-bookworm` | `buildpack-deps:bookworm` (Debian 12) | Python 3.12, pip, gcc/g++, curl, wget, gnupg, ca-certificates, build headers | ~1.0 GB |
| `python:3.12-slim-bookworm` | `debian:bookworm-slim` | ca-certificates, netbase, tzdata only | ~130 MB |
| `python:3.12-alpine` | Alpine Linux | Python 3.12, pip; glibc unavailable (musl) | ~60 MB |
| `python:3.12-windowsservercore` | Windows | N/A — wrong OS | N/A |

**Why `python:3.12-bookworm` wins over `slim`:** The slim variant excludes `build-essential`, `curl`, `wget`, `gnupg`, and dev headers, so we would have to reinstall them all for Octave to compile and for SSH setup. The full image bundles them via `buildpack-deps:bookworm`, making the installed apt list shorter.

**Why `bookworm` over `alpine`:** Octave and `octave-statistics` are not available in Alpine's apk repository. Alpine uses musl libc instead of glibc, which breaks Octave's Fortran-linked BLAS/LAPACK. Not viable.

**Debian 12 (Bookworm) compatibility:** Ubuntu 24.04 Noble Numbat is also Debian-based with very similar apt package availability, so the existing `octave` and `octave-statistics` apt installs work identically.

---

## Packages Pre-installed in `python:3.12-bookworm` (can be removed from the manifest)

Via `buildpack-deps:bookworm`:

| Package | Status |
|---|---|
| `python3`, `python-is-python3` | ✅ Pre-installed as `/usr/local/bin/python` |
| `python3-pip` | ✅ Pre-installed as `/usr/local/bin/pip` |
| `python3-dev` | ✅ C headers bundled via buildpack-deps |
| `build-essential` (gcc, g++, make) | ✅ Pre-installed via buildpack-deps |
| `curl` | ✅ Pre-installed via buildpack-deps |
| `wget` | ✅ Pre-installed via buildpack-deps |
| `gnupg` | ✅ Pre-installed via buildpack-deps |
| `ca-certificates` | ✅ Pre-installed explicitly |

---

## Packages Still Needed via `apt-get`

| Package | Why still needed |
|---|---|
| `openssh-server` | SSH daemon for port-forward access |
| `rsync` | Code sync from Mac dev machine |
| `tmux` | Session multiplexing in the pod |
| `htop` | Process monitoring |
| `octave` | Shambhala2 CuBlock algorithm |
| `octave-statistics` | Octave stats package dependency |

`awscli` is currently installed via `apt-get`. Since `boto3` and `botocore` are already in `requirements.txt` (covering all S3 operations), `awscli` is used only for the `aws` CLI binary. It should be **removed from apt-get and added to `requirements.txt`** as `awscli>=2.13` — this ensures a current version (apt's `awscli` on Bookworm is pinned to an old 1.x release).

---

## Required Changes to `shambhala-pod.yaml`

### Step 1 — Change the base image (line 43)

```yaml
# Before:
image: ubuntu:24.04

# After:
image: python:3.12-bookworm
```

### Step 2 — Trim the system package install (startup script, section 2)

Remove from the `apt-get install` call in section 2:
- `build-essential`
- `python3`
- `python3-pip`
- `python3-dev`
- `python-is-python3`
- `curl` (already present)
- `wget` (already present)
- `gnupg` (already present)
- `ca-certificates` (already present)
- `awscli` (moving to pip)

Reduced apt-get install line:
```bash
apt-get install -y -qq openssh-server rsync tmux htop octave octave-statistics
```

> Note: `apt-get update -qq` in section 1 is still required to refresh Bookworm's index before installing octave.

### Step 3 — Drop `--break-system-packages` from pip install (section 3)

The official Python image's pip is not a system-managed package, so the PEP 668 restriction does not apply. Remove the flag:

```bash
# Before:
pip3 install --break-system-packages --no-cache-dir -r /etc/shambhala-deps/requirements.txt

# After:
pip install --no-cache-dir -r /etc/shambhala-deps/requirements.txt
```

### Step 4 — Update `requirements.txt` ConfigMap to add awscli

In the `shambhala-deps` ConfigMap, add:
```
awscli>=2.13
```

### Step 5 — Update verification commands (section 3, optional but clean)

```bash
# Before:
python3 --version && pip3 --version

# After:
python --version && pip --version
```

`python3` still works (symlinked), but `python` is the canonical binary in the official image.

### Step 6 — Remove the `octave --version` check (optional cleanup)

The version print currently appears in section 2. It can stay for diagnostic purposes or be folded into section 3's verification step — no correctness impact either way.

---

## Summary of Changes

| Location | Change |
|---|---|
| `image:` field | `ubuntu:24.04` → `python:3.12-bookworm` |
| Section 1 `apt-get install` | Remove `build-essential python3 python3-pip python3-dev python-is-python3 curl wget gnupg ca-certificates awscli` |
| Section 3 `pip install` | Remove `--break-system-packages`; use `pip` not `pip3` |
| `requirements.txt` ConfigMap | Add `awscli>=2.13` |
| Version check line | `python3`→`python`, `pip3`→`pip` (cosmetic) |

---

## Expected Impact

- **Faster startup:** Section 2 apt install goes from ~10 packages to 6. The Python toolchain (previously the slowest part to configure on bare Ubuntu) is zero-cost.
- **No functional changes:** Same Octave version, same Python packages, same SSH/rsync setup.
- **Reproducibility:** Pinning to `python:3.12-bookworm` gives a stable minor-version image; the tag does not float to 3.13+.

---

## To-Do List

### 1. Edit `harmonization_scripts/k8s/shambhala-pod.yaml`

- [x] **1.1** Change `image: ubuntu:24.04` → `image: python:3.12-bookworm` (line 43).
- [x] **1.2** In the section-1 apt block (SSH setup), remove `curl wget gnupg ca-certificates` from the `apt-get install` line — they are pre-installed. Keep `openssh-server rsync tmux htop`.
  - Current line: `apt-get install -y -qq openssh-server rsync curl wget gnupg ca-certificates tmux htop`
  - New line: `apt-get install -y -qq openssh-server rsync tmux htop`
- [x] **1.3** In the section-2 apt block (system libs), replace the entire `apt-get install` line:
  - Remove: `build-essential python3 python3-pip python3-dev python-is-python3 awscli`
  - Remove: `curl wget gnupg ca-certificates` (already cleaned in 1.2 above; ensure they are not duplicated here)
  - Keep: `octave octave-statistics`
  - New line: `apt-get install -y -qq octave octave-statistics`
- [x] **1.4** Update the version-check line in section 2:
  - `python3 --version && pip3 --version` → `python --version && pip --version`
- [x] **1.5** In the section-3 pip install line, remove `--break-system-packages`:
  - `pip3 install --break-system-packages --no-cache-dir -r /etc/shambhala-deps/requirements.txt` → `pip install --no-cache-dir -r /etc/shambhala-deps/requirements.txt`
- [x] **1.6** Update the Python import verification line in section 3 to use `python` instead of `python3`:
  - `python3 -c "import boto3, ..."` → `python -c "import boto3, ..."`

### 2. Edit the `shambhala-deps` ConfigMap in the same file

- [x] **2.1** Add `awscli>=2.13` to the `requirements.txt` data block (after `psutil>=5.9`).
  - Reason: `awscli` is moving from apt (pinned to old 1.x on Bookworm) to pip to get a current 2.x release.

### 3. Delete the old pod and re-apply the manifest

- [ ] **3.1** Delete the running pod: `kubectl delete pod shambhala-bench -n <your-namespace>`
- [ ] **3.2** Apply the updated manifest: `kubectl apply -f harmonization_scripts/k8s/shambhala-pod.yaml`
- [ ] **3.3** Watch pod startup logs to confirm all sections complete without errors:
  `kubectl logs -f shambhala-bench -n <your-namespace>`

### 4. Verify the pod is functional

- [ ] **4.1** Port-forward SSH: `kubectl port-forward pod/shambhala-bench 2222:22 -n <your-namespace>`
- [ ] **4.2** SSH into the pod: `ssh -p 2222 root@localhost`
- [ ] **4.3** Inside the pod, verify the Python version: `python --version` — expect `Python 3.12.x`.
- [ ] **4.4** Verify pip and all required packages: `python -c "import boto3, pandas, numpy, sklearn, psutil; print('OK')"`.
- [ ] **4.5** Verify Octave is functional: `octave --no-gui --eval "disp('octave ok')"`.
- [ ] **4.6** Verify the `matlab` wrapper script still exists and strips unsupported flags: `matlab -nosplash -nodesktop --eval "disp('wrapper ok')"`.
- [ ] **4.7** Verify the `aws` CLI binary is available (from pip-installed awscli): `aws --version` — expect `aws-cli/2.x.x`.
- [ ] **4.8** Rsync a small test file from the Mac dev machine to confirm the sync workflow is unaffected.

### 5. Smoke-test the Shambhala pipeline

- [ ] **5.1** Sync the `Shambhala_containerized/` directory into `/app/Shambhala_containerized` on the pod.
- [ ] **5.2** Run the toy fixture test to confirm end-to-end Octave bridge still works:
  ```bash
  cd /app/Shambhala_containerized
  python run_shambhala.py \
      --input tests/fixtures/input_10samples.csv \
      --P tests/fixtures/P0_small.csv \
      --Q tests/fixtures/Q0_small.csv \
      --output /tmp/harmonized_test.csv \
      --octave-bin /usr/bin/octave \
      --n-workers 1 --random-seed 42 --log-level DEBUG
  ```
- [ ] **5.3** Confirm output matches reference within tolerance (`atol=1e-4`).

### 6. Update documentation

- [x] **6.1** Update `CLAUDE.md` in `Shambhala_containerized/` — changed `ubuntu:24.04` to `python:3.12-bookworm` in the K8s pod table entry.
- [x] **6.2** Mark this plan's to-do items as complete once verified.
