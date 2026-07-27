# Pod Python/pip Missing — Root Cause Analysis and Fix Plan

**Date:** 2026-05-15  
**Symptom:** Inside `shambhala-bench`: `python` and `pip3` commands not found despite
the startup log printing `[startup] Python packages installed`.

---

## 1. What the Logs Show

Relevant excerpt from `logs/shambhala_pod_logs_260515.txt`:

```
[startup] SSH ready - port-forward and connect now to watch progress   ← step 1 succeeded
E: Unable to locate package liboctave-dev                              ← step 2 FAILED
[startup] System libs + Octave installed                               ← false positive
bash: line 24: pip3: command not found                                 ← step 3 FAILED
[startup] Python packages installed                                    ← false positive
[startup] Environment ready. Sync code with rsync and run.            ← false positive
```

Steps 2 and 3 both failed silently. The `[startup]` echo lines ran regardless.

---

## 2. Root Cause #1 — `liboctave-dev` Does Not Exist in Ubuntu 24.04

The step 2 command in `shambhala-pod.yaml` is:

```bash
apt-get install -y -qq build-essential python3 python3-pip python3-dev python-is-python3 \
    octave octave-statistics liboctave-dev
```

`liboctave-dev` existed in Ubuntu 22.04 but was **removed from Ubuntu 24.04**. When `apt-get
install` encounters a package it cannot resolve, it aborts the entire command with exit code 100
and installs **nothing** from that command (apt resolves the full dependency graph before writing
anything to disk).

This means none of the following were installed:
- `python3-pip` → explains `pip3: command not found`
- `python-is-python3` → explains `python: command not found` (the symlink `/usr/bin/python → python3` is never created)
- `octave` → Octave binary not present
- `octave-statistics` → Octave statistics package not present
- `build-essential`, `python3-dev` → also absent

`python3` itself IS present on the pod because it was pulled in as a transitive dependency of
packages installed in step 1 (`openssh-server`, `tmux`, etc. depend on python3). This is why
`python3` works but `python` and `pip3` do not.

**`liboctave-dev` was never needed.** The Shambhala pipeline calls Octave as a subprocess via
stdin/stdout pipes (`subprocess.Popen`). It never links against liboctave C++ headers. The
package was included by mistake.

---

## 3. Root Cause #2 — No Error Checking in Startup Script

The pod startup script (`args` field in the YAML) does not begin with `set -e`. Each step's echo
line runs unconditionally whether or not the previous command succeeded. This produces false
`[startup] ... installed` messages even when the install failed, masking the failure completely.

The symptom: the pod reaches `[startup] Environment ready.` and stays Running in K8s, but is
non-functional.

---

## 4. Root Cause #3 — pip Install Target is Missing

Step 3 runs:
```bash
pip3 install --break-system-packages --no-cache-dir -r /etc/shambhala-deps/requirements.txt
```

Even if `liboctave-dev` is removed and `python3-pip` installs correctly, the
`--break-system-packages` flag is required on Ubuntu 24.04 because Python 3.12 marks the system
site-packages as externally managed. This flag is already present in the YAML — no change needed
here once the apt-get succeeds.

---

## 5. Implementation Plan

### Change 1 — Remove `liboctave-dev`; add `set -e`

Edit `harmonization_scripts/k8s/shambhala-pod.yaml`, pod `args` field.

**Add `set -e` as the very first line** of the bash script so any future failure aborts startup
immediately with a visible error instead of silently continuing.

**Remove `liboctave-dev`** from the step 2 apt-get command.

Before:
```bash
apt-get install -y -qq build-essential python3 python3-pip python3-dev python-is-python3 \
    octave octave-statistics liboctave-dev
```

After:
```bash
apt-get install -y -qq build-essential python3 python3-pip python3-dev python-is-python3 \
    octave octave-statistics
```

### Change 2 — Add Octave verification step

After installing Octave, add a line that verifies the binary exists and fails loudly if not:
```bash
octave --version || { echo "[startup] ERROR: octave binary missing"; exit 1; }
```

### Change 3 — Add python verification step

After pip install, verify the critical binaries:
```bash
python3 --version && pip3 --version || { echo "[startup] ERROR: python/pip missing"; exit 1; }
which python   || { echo "[startup] WARNING: python symlink missing"; }
```

---

## 6. Updated Pod Args (Full Step 2 Block)

Replace the current step 2 block with:

```bash
# -- 2. System libs for Python and Octave (no R, no liboctave-dev) --
apt-get install -y -qq build-essential python3 python3-pip python3-dev python-is-python3 \
    octave octave-statistics
octave --version
python3 --version && pip3 --version
printf '#!/bin/bash\nargs=()\nfor a in "$@"; do\n  case "$a" in\n    -nodesktop|-nosplash|-nodisplay) ;;\n    *) args+=("$a") ;;\n  esac\ndone\nexec octave --no-gui "${args[@]}"\n' \
    > /usr/local/bin/matlab && chmod +x /usr/local/bin/matlab
echo "[startup] System libs + Octave installed"
```

The only substantive change is removing `liboctave-dev`. Everything else is the same or adds
verification.

---

## 7. Why No Other Changes Are Needed

| Concern | Status |
|---|---|
| `python3-pip` provides `pip3` | Correct — no change needed |
| `python-is-python3` provides the `python` symlink | Correct — no change needed |
| `--break-system-packages` in pip install | Already present — correct for Ubuntu 24.04 |
| Octave statistics package (`octave-statistics`) | Still present and correct |
| `octave` binary path inside pod | `/usr/bin/octave` — correct, no change needed |
| `matlab` shim script | Unchanged — correct |
| AWS credentials, PVC, resource limits | Unchanged — correct |

---

## 8. Validation Checklist (After Redeploying)

After applying the fix and running `kubectl apply -f ... && kubectl logs -f ...`:

```
# Expected in logs:
[startup] SSH ready - port-forward and connect now to watch progress
octave, version X.Y.Z
Python 3.12.x
pip 24.x ...
[startup] System libs + Octave installed
[startup] Python packages installed
[startup] Environment ready. Sync code with rsync and run.
```

After SSH into pod:
```bash
python3 --version        # Python 3.12.x
python  --version        # Python 3.12.x  (via python-is-python3 symlink)
pip3    --version        # pip 24.x
octave  --version        # GNU Octave, version x.y.z
which   matlab           # /usr/local/bin/matlab
```

Then proceed with:
```bash
rsync ... /app/Shambhala_containerized/
cd /app/Shambhala_containerized && pip3 install --break-system-packages -e .
```
