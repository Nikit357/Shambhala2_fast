# K8s Pod Operations — Shambhala Benchmark

This guide covers deploying and operating the `shambhala-bench` Kubernetes pod.
The pod runs the Shambhala cross-product benchmark (54 jobs → 1,512 S3 outputs) with no R dependency.

**Startup time:** ~5–10 min (no R packages to install).

---

## Prerequisites

- `kubectl` configured for your Kubernetes cluster.
- `~/.aws/credentials` present with S3 read/write access to your bucket.
- PVC `shambhala-pvc` already exists in your namespace.
- `SHAMBHALA_S3_BUCKET` exported — see `.env.example` in the repository root.

---

## One-time setup

Create the AWS credentials secret (only needed once per cluster):

```bash
kubectl create secret generic aws-credentials \
    --from-file=credentials=$HOME/.aws/credentials \
    --from-file=config=$HOME/.aws/config \
    --namespace <your-namespace>
```

---

## Deploy the pod

```bash
kubectl apply -f harmonization_scripts/k8s/shambhala-pod.yaml -n <your-namespace>
```

Check that the pod was created and is starting:

```bash
kubectl get pod shambhala-bench -n <your-namespace>
```

Monitor startup logs (wait for `[startup] Environment ready.`):

```bash
kubectl logs -f shambhala-bench -n <your-namespace>
```

Startup completes in ~5–10 min. SSH becomes available within ~1 min (before Octave finishes installing).

---

## SSH access

Port-forward and connect:

```bash
kubectl port-forward pod/shambhala-bench 2223:22 -n <your-namespace> &
ssh -p 2223 -i ~/.ssh/id_ed25519 -o StrictHostKeyChecking=no root@localhost
```

Or add to `~/.ssh/config` for convenience:

```
Host shambhala-pod
    HostName localhost
    Port 2223
    User root
    IdentityFile ~/.ssh/id_ed25519
    StrictHostKeyChecking no
```

Then connect with: `ssh shambhala-pod`

---

## Sync code from Mac to pod

```bash
rsync -avz --progress \
    -e "ssh -p 2223 -i ~/.ssh/id_ed25519 -o StrictHostKeyChecking=no" \
    /Users/user890/Desktop/fl_subset/shambhala_adoption/Shambhala_containerized/ \
    root@localhost:/app/Shambhala_containerized/
```

---

## Install the Shambhala package inside the pod

After syncing code:

```bash
cd /app/Shambhala_containerized
pip install -e .
```

---

## Choosing `--n-workers` and `--n-shambhala-workers`

These two parameters control parallelism at two independent levels, and must be set correctly to avoid Octave timeout failures.

### Two-level timeout hierarchy

The dispatcher operates with two separate timeout limits:

| Argument | Controls | Default |
|---|---|---|
| `--timeout-s` | How long the dispatcher waits for a full `run_shambhala_job.py` subprocess | 21600 s (6 h) |
| `--shambhala-timeout-s` | How long each individual Octave call inside the job may run | inherits `--timeout-s` if not set |

**Both must be set** when the batch size (NH) is large. If `--shambhala-timeout-s` is omitted it falls back to `--timeout-s`, which is the correct default after the 2026-05-17 fix. Before this fix, the inner timeout was hardcoded to 7200 s regardless of `--timeout-s`.

### How NH (samples per Octave call) is determined

```
NH = n_S0_samples / --n-shambhala-workers
```

S0 strict has 7174 samples. With `--n-shambhala-workers 5`, NH = 1435. With 16, NH ≈ 449.

CuBlock runs 30 iterations of k-means + cubic fitting on a matrix of `n_genes × (NH + NP)` in interpreted (non-JIT) Octave. Computation time scales roughly linearly with NH + NP. NH=1435 exceeded 2 hours per batch; NH≈449 is estimated at 30–60 minutes.

### Total simultaneous Octave processes

```
total_octave_processes = --n-workers × --n-shambhala-workers
```

This must not exceed the pod vCPU count. The shambhala pod has **16 vCPU**. Each Octave process is single-threaded.

### Recommended configurations for the 16-vCPU shambhala pod

| `--n-workers` | `--n-shambhala-workers` | Total Octave processes | NH (strict S0) | Notes |
|---|---|---|---|---|
| 1 | 16 | 16 | 449 | **Recommended start** — saturates CPU, no contention between jobs |
| 2 | 8 | 16 | 898 | 2 jobs in parallel; NH larger so each batch takes longer |
| 1 | 30 | 30 | 239 | Over-subscribes CPUs slightly; may be faster overall if NH=449 is bottlenecked by memory |
| 5 | 5 | 25 | 1435 | **Do not use** — NH too large, all batches timed out in the 2026-05-17 run |

### Observed timing data

| `--n-shambhala-workers` | NH (strict S0, 7174 samples) | Elapsed per Octave batch | Date tested |
|---|---|---|---|
| 5 | 1435 | >7200 s (timed out) | 2026-05-17 |
| 16 | 449 | TBD | — |

Fill in the TBD entry after the first successful test job with `--n-shambhala-workers 16`.

---

## Run the benchmark

Recommended command for the 16-vCPU pod (run from `/app/Shambhala_containerized/` inside the pod):

```bash
nohup python harmonization_scripts/run_shambhala_parallel.py \
    --n-workers 1 \
    --n-shambhala-workers 16 \
    --octave-bin /usr/bin/octave \
    --skip-if-exists \
    --memory-limit-gb 10.0 \
    --timeout-s 7200 \
    --random-seed 42 \
    > /workspace/shambhala_run.log 2>&1 &
```

`--shambhala-timeout-s` is omitted here; it inherits `--timeout-s 7200` automatically.

To set the inner Octave timeout independently (e.g. if you want a longer outer safety net):

```bash
nohup python harmonization_scripts/run_shambhala_parallel.py \
    --n-workers 1 \
    --n-shambhala-workers 16 \
    --octave-bin /usr/bin/octave \
    --skip-if-exists \
    --memory-limit-gb 10.0 \
    --timeout-s 14400 \
    --shambhala-timeout-s 7200 \
    --random-seed 42 \
    > /workspace/shambhala_run.log 2>&1 &
```

The last example by Daniil:

```bash
python harmonization_scripts/run_shambhala_parallel.py  --n-workers 1  --n-shambhala-workers 30  --octave-bin /usr/bin/octave   --skip-if-exists  --memory-limit-gb 10.0    --timeout-s 216000  --shambhala-timeout-s 72000   --random-seed 42  --methods shambhala_P0std_Q0std,shambhala_NBGPL570_Q0std,shambhala_NBKass_Q0std,shambhala_NBRNAseq_Q0std   > shambhala_run_260518.log
```


To validate a single job before the full run:

```bash
python harmonization_scripts/run_shambhala_job.py \
    --imp strict \
    --method shambhala_P0std_Q0std \
    --n-shambhala-workers 16 \
    --octave-bin /usr/bin/octave \
    --shambhala-timeout-s 7200 \
    --random-seed 42 \
    --out-json /tmp/test_job.json
```

Watch for `Octave batch done: NH=..., NP=..., elapsed=...s.` lines in the output to confirm per-batch timing.

---

## Check progress

```bash
# Live log
tail -f /workspace/shambhala_run.log

# Count S3 outputs (target: 1296)
aws s3 ls s3://your-bucket/FL_batch_correction/exp/ \
    | grep shambhala | wc -l

# CPU / RAM usage
htop
```

---

## Teardown

```bash
kubectl delete pod shambhala-bench -n <your-namespace>
```

The PVC `/workspace` and all S3 outputs persist after pod deletion.

---

## Speed-Up Flags

All speed-up flags are available in the pod. Full documentation in
`harmonization_scripts/README.md` (Performance Tuning section).

**Key flags for NBGPL570 (NP=250, the slowest variant):**

```bash
# Recommended for production: Python QN reference (2–4× speedup, < 0.05 max deviation)
--precompute-qn-reference

# Maximum speed: Python QN + synthetic P + Python CuBlock (~10–20× speedup)
--precompute-qn-reference --synthetic-cublock-p --python-cublock

# Validate a new combination before full run:
pytest tests/speed_up_tests/ -v -s 2>&1 | tee /workspace/speed_up_test_results.txt
```

**Dependencies required for speed-up flags** (already in ConfigMap pip deps):
- `qnorm>=0.4.1` — used by `--precompute-qn-reference` and `--python-cublock`
- `scikit-learn>=1.3` — used by `--precompute-cublock-clusters` and `--python-cublock`
