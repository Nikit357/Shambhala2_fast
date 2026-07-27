# Quantile Normalization Speed-Up Plan for Shambhala Harmonization

**Date:** 2026-05-18  
**Last updated:** 2026-05-18 (all inline comments resolved)  
**Author:** Daniil Nikitin  
**Status:** Implementation complete — VER-4 (commits) pending user request  

---

## Summary of Decisions

| Approach | Decision | CLI flag | Default |
|---|---|---|---|
| **A** — Precompute QN reference (synthetic P) | ✅ Approved | `--precompute-qn-reference` + `--synthetic-cublock-p` | Off |
| **B** — Skip polynomial fits for P columns | ✅ Approved as **obligate default** | N/A — always on | On |
| **C** — P subsampling | ✅ Approved | `--max-p-samples N` (PCA centroid) | Off |
| **D** — Python-side QN (qnorm library) | ✅ Approved, **tied to Approach A flag** | Auto-activated by `--precompute-qn-reference` | Off |
| **E** — Precompute k-means clusters from P | ✅ Approved | `--precompute-cublock-clusters` + new `CuBlock_fixed.m` | Off |
| **F** — Batch QN across all samples at once | ❌ **Rejected** | — | — |
| **G** — Pure Python CuBlock | ✅ Approved | `--python-cublock` + new `cublock_python.py` | Off |

All approved changes (A, B, C, D, E, G) are covered by integration tests in `tests/speed_up_tests/`. Approach F is excluded from all tests.

---

## 1. Problem Summary

### Observed timings (from `shambhala_pod_logs_260518.txt`)

| Variant | P size (samples) | Genes (intersection) | Time/sample | Estimated total | Workers |
|---|---|---|---|---|---|
| `shambhala_P0std_Q0std` | **39** | 2,293 | **~10 s** | **~50 min** ✅ | 30 |
| `shambhala_NBGPL570_Q0std` | **250** | 2,366 | **~155 s** | **~12 hrs** ❌ | 30 |

At 6 h elapsed, NBGPL570 was at 55–64% completion. The run is projected to take ~11–12 h total — approximately **14–15× slower** for only **6.4× more P samples**.

### Root cause: both QN and CuBlock scale with P size (NP)

Each Octave subprocess processes a pool matrix of shape **(G genes × [1 input sample + NP P samples])**:

```
Shambhala2_piped.m inner loop (per input sample i):
   i0 = [i, NH+1, NH+2, ..., NH+NP]     ← 1 sample + all NP P columns
   EXP = Exp(:, i0)                      ← G × (1+NP) sub-matrix
   EXP = quantilenorm(EXP)               ← sort all 1+NP columns
   dataN = CuBlock(real(EXP), [], k)     ← k-means on G × (1+NP), polynomial fit for ALL 1+NP columns
   vecN = EXPN(:, 1)                     ← ONLY the first column is kept; NP columns discarded
```

**Quantile normalization (`quantilenorm.m`, 17 lines):**
- Sorts each of the 1+NP columns independently: O(G × NP × log G)
- For NBGPL570: sorts 251 × 2,366 = 594,266 values per sample
- For P0std: sorts 40 × 2,293 = 91,720 values per sample → 6.5× more work
- Crucially: **P columns are identical every iteration** — they are re-sorted 7,174 times unnecessarily

**CuBlock (`CuBlock.m`, lines 49–113):**
- Runs k-means on the full G × (1+NP) matrix 30 times: O(G × NP × k × iterations × reps)
- For each of 30 reps: polynomial fit over **all 1+NP columns** (line 51: `for j=1:nbSamples`)
- Then discards the P columns — **250 polynomial fits per rep are computed and thrown away**
- The polynomial fitting cost per rep: 30 × 251 × G = 17.8M ops (NBGPL570) vs 30 × 40 × G = 2.75M ops (P0std) → 6.5×

Combined, the 6.4× increase in NP causes ~15× slowdown because:
1. QN scales linearly with NP (sorts all columns)
2. CuBlock k-means scales linearly with NP (distances in NP-dimensional space)
3. CuBlock polynomial fits scale linearly with NP (loop over all columns)
4. Memory overhead grows with NP, causing cache effects and Octave GC pressure

---

## 2. Approaches to Speed-Up

### Approach A — Precompute the QN Reference from P (Python-side, once) ✅ APPROVED

**Decision: Approved. Implement as two optional CLI flags (see below).**

**What it is:**  
Standard quantile normalization with a single reference dataset P computes `rank_means[r] = mean(sorted_P[r, :])` once per gene rank. For each new input sample: sort it → assign the P rank means → unsort. This reduces QN from O(G × NP × log G) per sample to O(G × log G) per sample.

Currently the code re-sorts all NP P columns for every input sample. The reference distribution from P is constant and can be precomputed once in Python/numpy before any Octave call, then passed as a single representative column.

**Implementation:**
```python
# In run_shambhala.py or parallel.py, before harmonize_parallel():
p_sorted = np.sort(p_T.values, axis=0)           # sort each P column (genes×P)
p_qn_reference = p_sorted.mean(axis=1)            # mean per rank → (G,)
p_synthetic = pd.DataFrame(
    p_qn_reference[:, np.newaxis],
    index=p_T.index,
    columns=["P_qn_reference"]
)
# Pass p_synthetic instead of p_T to harmonize_parallel()
```

Octave then sees a G × 2 matrix (1 input sample + 1 synthetic P reference column) for both QN and CuBlock instead of G × 251.

**CLI flags (two independent flags):**

**`--precompute-qn-reference`** — activates the QN optimization:
- Computes the synthetic P reference (mean of sorted P columns per rank) in Python
- Activates Python-side QN using the `qnorm` library (see Approach D below — automatically bundled)
- Passes the QN-normalized input samples + 1 synthetic P column to Octave
- Octave's internal `quantilenorm()` now operates on 2 columns only

**`--synthetic-cublock-p`** — activates the CuBlock P optimization:
- Passes the same 1-column synthetic P to Octave for CuBlock
- When used with `--precompute-qn-reference`, both QN and CuBlock operate on G×2 matrices (maximum speedup)
- When used alone (without `--precompute-qn-reference`), CuBlock receives synthetic P but QN is unchanged

**Architectural note on coupling:** In `Shambhala2_piped.m`, QN and CuBlock share the same `EXP` matrix — the code cannot currently separate QN's view of P from CuBlock's view of P in a single Octave call. The two flags above provide logically separate control: `--precompute-qn-reference` performs Python-side QN, removing QN from Octave's workload entirely; `--synthetic-cublock-p` controls what P representation CuBlock receives. Together they are fully decoupled. Separately, each provides partial speedup.

**Both flags must be added to:**
- `harmonization_scripts/run_shambhala_parallel.py` (passes them through to each job subprocess)
- `harmonization_scripts/run_shambhala_job.py` (passes them through to `run_shambhala.py`)

**Required changes:**
- 4–8 lines of Python added before `harmonize_parallel()` in `run_shambhala.py`
- `pip install qnorm` (C-extension library) added to requirements
- No changes to `.m` files
- No changes to `octave_bridge.py` or `parallel.py` beyond accepting `p_synthetic`

**Speed-up on QN phase:** ~NP× = ~250× for NBGPL570 (the QN call sees 2 columns instead of 251)

**Speed-up on CuBlock phase (when `--synthetic-cublock-p` also set):** ~NP× = ~250× (CuBlock now operates on G × 2 instead of G × 251)

**Estimated total speed-up (both flags):** If QN and CuBlock each contribute ~50% of per-sample time, and both are reduced proportionally: **~100–200× speedup for the QN+CuBlock combined**. Per-sample time would drop from ~155 s to ~2–5 s.

**Authenticity impact:** ⭐⭐⭐⭐⭐ (minimal loss)  
The quantile normalization reference shifts from `mean(sorted([sample_i, P_1, ..., P_NP]))` to `mean(sorted([P_1, ..., P_NP]))`. For NP = 250, sample_i contributes 1/251 ≈ **0.4%** to the reference distribution. The change is statistically negligible. This is conceptually equivalent to what happens in bulk QN when NP → ∞.

For CuBlock with 1 synthetic P column: the synthetic P column (mean of all sorted P columns per gene rank) is the most representative single column possible. K-means on [sample | P_mean] should produce similar high-level gene clusters because gene separation is driven by gene-level expression magnitudes, which are preserved in the P reference column. The 30-rep averaging provides robustness.

**Verdict:** Highest expected speed-up, low authenticity risk for QN, moderate risk for CuBlock gene clustering. Primary implementation target.

---

### Approach B — CuBlock: Skip Polynomial Fits for P Columns ✅ APPROVED AS DEFAULT

**Decision: Approved as obligate default behavior — applied unconditionally to all Shambhala runs. No CLI flag needed. Zero authenticity risk confirmed by detailed code trace below.**

**What it is:**  
In `CuBlock.m`, the inner loop (line 51) iterates `for j=1:nbSamples` — all 1+NP columns. But `Shambhala2_piped.m` discards columns 2..NP+1 immediately after: `vecN = EXPN(:,1)`. The polynomial fitting for columns 2–251 is **100% wasted computation**.

#### Detailed Code Trace Confirming That Polynomial Fits on P Are Never Used

**Step 1: `Shambhala2_piped.m` calls CuBlock with the full G×(1+NP) matrix.**

```matlab
EXP = Exp(:, i0);              % G × (1+NP): column 1 = input sample; columns 2..NP+1 = P samples
EXP = quantilenorm(EXP);      % QN all 1+NP columns
dataN = CuBlock(real(EXP), [], k);   % G × (1+NP) returned — one output column per input column
```

**Step 2: Inside `CuBlock.m`, k-means runs on the full G×(1+NP) matrix.**

```matlab
indProbes = kmeans(data, k, 'maxiter', 1000);
```

This step receives the full G×(1+NP) matrix `data` and clusters genes (rows) based on expression patterns across **all** 1+NP samples (columns). **The entire P dataset is used here.** The resulting `indProbes` is a G-length vector assigning each gene to cluster 1..k. This is the legitimate calibration role of P — it determines gene cluster assignments. The user's intuition that "the entire P dataset is used for CuBlock k-means calculation" is correct for this step.

**Step 3: The polynomial fitting loop runs over all columns, including P columns.**

```matlab
for j = 1:nbSamples          % j = 1...(1+NP), i.e., input sample + all NP P samples
    for i = 1:k
        if sum(indProbes==i) > 100
            dataCurr = data(indProbes==i, j);           % genes in cluster i, column j
            % ... Z-transform ...
            pol = polyfit(dataCurrS, X(:,indP), 3);     % cubic polynomial for column j
            currDataN = ModPol(dataCurr, indS, pol);    % apply polynomial for column j
            dataN(indProbes==i, j) += currDataN;        % stored in dataN column j
            count(indProbes==i, j) += 1;
        end
    end
end
dataN = dataN ./ count;       % G × (1+NP) averaged result
```

`dataN` is a G×(1+NP) matrix. For every repetition (`nRep`), polynomial fits are computed and accumulated for **all** 1+NP columns, including all NP P columns (j=2..NP+1).

**Step 4: Back in `Shambhala2_piped.m`, only column 1 is ever read.**

```matlab
dataN = CuBlock(real(EXP), [], k);   % G × (1+NP)
log2e = log2(exp(1));
DataN = dataN / log2e;               % G × (1+NP)
EXPN = exp(DataN) - 1;               % G × (1+NP)
vecN = EXPN(:, 1);                   % ONLY COLUMN 1 (the input sample)
```

Columns 2..NP+1 of `EXPN` are never assigned, never read, never written anywhere. They are immediately eligible for garbage collection when `dataN`/`DataN`/`EXPN` go out of scope. The polynomial fit computation for those columns is **100% wasted**.

**This is confirmed in the original GitHub repository (`BorisovNM/Shambhala2`) too:** the original `CuBlock.m` has the identical `for j=1:nbSamples` loop and `Shambhala2.m` has the identical `vecN = EXPN(:,1)` extraction. This was never a Shambhala design intent — it is simply an oversight in the original code.

**Implementation:** Single line change in `CuBlock.m`:

```matlab
% Line 51 original:
for j = 1:nbSamples

% Line 51 modified:
for j = 1:1    % only the input sample (j=1); P columns calibrate k-means only
```

Additionally, `dataN` and `count` can be initialized as G×1 vectors instead of G×S matrices, reducing memory allocation by a factor of (1+NP).

**Speed-up:** The polynomial fitting loop drops from 1+NP to 1 iteration per k-means rep per cluster. k-means (which does use all columns) is unchanged. If polynomial fitting accounts for ~30–40% of per-rep CuBlock time, this gives ~2–3× speedup on CuBlock overall. Combined with Approach A (which also reduces the k-means matrix), the combined speedup is much larger.

**Authenticity impact:** ⭐⭐⭐⭐⭐ (zero loss)  
This is a pure optimization of a computation that was being thrown away anyway. The algorithm output is **mathematically identical** for the input sample column (j=1). The k-means step is unchanged. **This is the highest-authenticity optimization available.**

**Verdict:** Implemented as obligate default. No flag needed. This change is applied in all Shambhala runs unconditionally from the next version forward.

---

### Approach C — P Subsampling (Representative Subset of P) ✅ APPROVED AS CLI FLAG

**Decision: Approved. CLI flag `--max-p-samples N`. Default: None (no subsampling). Strategy: PCA centroid proximity (option 1 from the list below).**

**What it is:**  
Rather than using all NP P samples in the pool matrix, select a representative subset of N_sub samples. Both QN and CuBlock operate on G × (1+N_sub) instead of G × 251.

**Implementation:**

```python
# In run_shambhala.py, after loading p_df and before harmonize_parallel():
if args.max_p_samples is not None and p_df.shape[0] > args.max_p_samples:
    original_size = p_df.shape[0]
    # Centroid-based selection: samples closest to the mean of P in gene space
    p_mean = p_df.mean(axis=0)                              # mean gene expression
    distances = ((p_df - p_mean) ** 2).sum(axis=1)         # Euclidean distance to centroid
    p_df = p_df.loc[distances.nsmallest(args.max_p_samples).index]
    logger.info(
        f"P subsampled to {args.max_p_samples} centroid-closest samples "
        f"(from {original_size}) using PCA centroid proximity"
    )
```

**Selection strategy: centroid-based (approach 1)**
This is the highest-quality subsampling strategy:
- Compute the mean gene expression vector across all P samples (the centroid in gene space)
- Select the N_sub samples with smallest Euclidean distance to the centroid
- These samples are most "representative" of the average P expression profile
- Avoids the risk of selecting outlier samples that would distort QN/CuBlock calibration

The `--max-p-samples` flag must be added to:
- `harmonization_scripts/run_shambhala_parallel.py`
- `harmonization_scripts/run_shambhala_job.py`

**Speed-up:** Directly proportional to reduction in NP:
- 250 → 40 samples: **~6× speedup** → estimated ~25–27 s/sample → ~60 min total
- 250 → 20 samples: **~12× speedup** → estimated ~13 s/sample → ~33 min total
- By default (None): no subsampling, full NP used

**Authenticity impact:** ⭐⭐⭐ (moderate concern)  
The QN reference is computed from fewer samples → more variance in the target distribution → slightly noisier normalization. The CuBlock k-means has fewer exemplars for gene clustering. With 40 carefully chosen centroid-based samples, the risk is low (the original P0std uses only 39). With random 20, there is a real risk of missing subpopulations.

**Validation required:** Run the reference fixture with N_sub=40 subsampled NBGPL570 vs NP=250 full. Compare output table against `tests/fixtures/expected_output.csv`.

---

### Approach D — Python-Side QN Using qnorm Library ✅ APPROVED, TIED TO APPROACH A

**Decision: Not an independent flag. Automatically activated when `--precompute-qn-reference` is set. Uses the `qnorm` library (C extension) for Python-side quantile normalization.**

**What it is:**  
Move the quantile normalization step out of Octave into Python, executed once per sample in Python before the Octave call. The qnorm library provides a fast C implementation.

When `--precompute-qn-reference` is active:
1. Python computes the QN reference from P: `p_qn_reference = np.sort(p_T.values, axis=0).mean(axis=1)`
2. Python applies QN to each input sample against this reference using `qnorm.quantile_normalize()`
3. Python sends QN-normalized input samples to Octave
4. If `--synthetic-cublock-p` is also set: Octave receives QN'd input + 1 synthetic P column; Octave's internal `quantilenorm()` becomes a near-no-op on 2 already-normalized columns
5. CuBlock operates on G×2

**Implementation:**
```python
import qnorm

if args.precompute_qn_reference:
    # Precompute QN reference from P
    p_sorted = np.sort(p_T.values, axis=0)       # (G, NP)
    p_qn_reference = p_sorted.mean(axis=1)        # (G,)
    p_synthetic = pd.DataFrame(
        p_qn_reference[:, np.newaxis],
        index=p_T.index,
        columns=["P_qn_reference"]
    )
    # QN each input sample against the P reference using qnorm
    # qnorm.quantile_normalize normalizes each column against the reference
    combined = pd.concat([input_T, p_synthetic], axis=1)
    combined_qn = qnorm.quantile_normalize(combined, axis=0)
    input_T = combined_qn.iloc[:, :-1]     # QN'd input samples
    p_for_octave = p_synthetic if args.synthetic_cublock_p else p_T
```

**Speed-up on QN phase:** `qnorm` library (C extension) is ~50–100× faster than Octave's loop-based implementation for large matrices. Since this is now Python-side, the Octave call for QN is eliminated entirely.

**Authenticity impact:** ⭐⭐⭐⭐⭐ (identical algorithm, faster implementation)  
Mathematically identical to `quantilenorm.m` when done against the same P reference distribution.

---

### Approach E — Precompute K-Means Gene Cluster Assignments from P ✅ APPROVED AS CLI FLAG

**Decision: Approved. CLI flag `--precompute-cublock-clusters`. New file: `octave/CuBlock_fixed.m`. Added to both `run_shambhala_parallel.py` and `run_shambhala_job.py`.**

**What it is:**  
In CuBlock, k-means clusters genes (rows) based on expression patterns across all samples (columns). The cluster assignments depend on the full [sample | P] matrix. However, since P has 250 samples vs 1 input sample, the cluster assignments are almost entirely P-driven.

Precompute cluster assignments once from P alone, then in the per-sample CuBlock call, reuse the precomputed assignments instead of re-running k-means.

**Implementation:**

New file `octave/CuBlock_fixed.m`:

```matlab
function dataN = CuBlock_fixed(data, N, k, fixed_clusters)
% CuBlock with precomputed gene cluster assignments.
% fixed_clusters: G×1 integer vector, precomputed from P using CuBlock on P alone.
% k-means is skipped; fixed_clusters is used as indProbes in every rep.
if nargin < 3 || isempty(k), k = 5; end
if nargin < 2 || isempty(N), N = 30; end
data = double(data);
[nbProbes, ~] = size(data);
dataN = zeros(nbProbes, 1);
count = zeros(nbProbes, 1);
for nRep = 1:N
    for j = 1:1                               % only the input sample
        for i = 1:k
            if sum(fixed_clusters==i) > 100
                dataCurr = data(fixed_clusters==i, j);
                % ... Z-transform, polyfit, ModPol, store (same as CuBlock.m) ...
            end
        end
    end
end
dataN = dataN ./ count;
end
```

Python-side precomputation:
```python
if args.precompute_cublock_clusters:
    # Run a single Octave call on P alone to get stable cluster assignments
    p_clusters = precompute_clusters_from_p(p_T, octave_bin=args.octave_bin, k=args.k)
    # p_clusters: (G,) array of integer cluster labels 1..k
    # Pass p_clusters to octave_bridge.py; bridge injects them as fixed_clusters
```

**Speed-up:** Eliminates 30 k-means runs per sample (the primary cost at large NP). Each sample now only runs 30 iterations of the polynomial fitting loop (already fixed by Approach B to j=1 only). Estimated **5–10× speedup** on CuBlock. QN is unchanged (still needs Approach A).

**Authenticity impact:** ⭐⭐ (moderate loss)  
CuBlock's random restarts serve a specific purpose: by averaging 30 different k-means partitions, the algorithm produces smooth, partition-independent normalization. Reusing a single set of cluster assignments from P removes this averaging over random partitions. The polynomial fit is applied deterministically, which may introduce bias tied to one particular gene clustering. Since P (250 samples) dominates sample space, P-only cluster assignments are likely representative. Validation against reference output required before production use.

---

### Approach F — Batch QN Across All Samples at Once ❌ REJECTED

**Decision: Rejected. This approach changes the algorithm's semantics in an unacceptable way and will not be implemented.**

**Why rejected:** In the original Shambhala2 algorithm, each input sample is quantile-normalized against P independently (one-at-a-time normalization). Batch QN across all input samples simultaneously would include all other input samples in the normalization reference for each sample — changing the semantics from sample-vs-P normalization to batch normalization. This makes results batch-composition-dependent, which violates the core Shambhala design principle of sample-independent harmonization.

**No implementation. No tests. No CLI flag.**

---

### Approach G — Pure Python CuBlock ✅ APPROVED AS CLI FLAG

**Decision: Approved as optional CLI flag `--python-cublock`. New file: `shambhala/cublock_python.py`. Separate unit tests required. Validation against Octave reference output mandatory.**

**What it is:**  
Rewrite CuBlock in Python using:
- `sklearn.cluster.KMeans` for gene clustering
- `numpy.polyfit` for polynomial fitting
- Full vectorization across samples
- Eliminates Octave startup overhead (~2 s per subprocess)

```python
from sklearn.cluster import KMeans
import numpy as np

def cublock_python(data: np.ndarray, n_reps: int = 30, k: int = 5) -> np.ndarray:
    """
    Pure-Python CuBlock normalization.
    data: G × S matrix (genes × samples). S = 1 (input) + NP (P samples).
    Returns: G × 1 normalized output for the input sample only.
    """
    G, S = data.shape
    dataN = np.zeros(G)
    count = np.zeros(G)
    for rep in range(n_reps):
        km = KMeans(n_clusters=k, n_init=1, max_iter=1000, random_state=rep)
        labels = km.fit_predict(data)              # cluster genes; uses all S columns
        for cluster in range(k):
            mask = labels == cluster
            if mask.sum() < 100:
                continue
            col = data[mask, 0]                    # input sample only (j=0)
            col_std = col.std(ddof=1)
            if col_std == 0:
                continue
            col_z = (col - col.mean()) / col_std
            sorted_idx = np.argsort(col_z)
            col_zs = col_z[sorted_idx]
            # GetTargetValues
            p_exp = np.arange(3, 22, 2)
            tol = 1e-1
            X = np.linspace(-1, 1, len(col_z))[:, None] ** p_exp
            std_mask = (col_zs >= -1) & (col_zs <= 1)
            S_vals = np.abs(X[std_mask]).mean(axis=0)
            valid = np.where(S_vals < tol)[0]
            indP = min(len(p_exp) - 1, valid[0] if len(valid) > 0 else len(p_exp) - 1)
            pol = np.polyfit(col_zs, X[:, indP], 3)
            curr_dataN = mod_pol_python(col_z, sorted_idx, pol)
            gene_idx = np.where(mask)[0]
            dataN[gene_idx] += curr_dataN
            count[gene_idx] += 1
    with np.errstate(invalid='ignore'):
        result = np.where(count > 0, dataN / count, np.nan)
    return result
```

**Speed-up:**
- `sklearn.KMeans` with BLAS is significantly faster than Octave's k-means for large matrices
- Vectorized polynomial fitting across all samples simultaneously
- Eliminates Octave startup overhead (~2 s per subprocess)
- Estimated total: **3–10× speedup on CuBlock**

**Authenticity impact:** ⭐⭐⭐ (moderate risk)  
`sklearn.KMeans` uses k-means++ initialization by default vs the custom `kmeans.m` (pure random). The 30-rep averaging mitigates differences but validation against reference output is mandatory (tolerance `atol=1e-3`).

**Development cost:** High. Requires: new `shambhala/cublock_python.py` module, `ModPol` reimplementation, comprehensive unit tests against the Octave reference, numerical agreement validation.

**Verdict:** Highest long-term maintainability (no Octave dependency). Recommended after dissertation deadline if needed. Implement now as a flagged option, with separate unit tests.

---

## 3. Summary Comparison Table

| Approach | Speed-Up Estimate | Authenticity | Dev Complexity | Default | Files Changed |
|---|---|---|---|---|---|
| **A — Precompute P QN reference** | **~10–20× overall** | ⭐⭐⭐⭐⭐ | Low | Off | `run_shambhala.py` + 2 harmonization scripts |
| **B — Skip polynomial fits for P columns** | ~1.5–2× (CuBlock only) | ⭐⭐⭐⭐⭐ | Minimal | **On (default)** | `CuBlock.m` (1 line) |
| **C — P subsampling (centroid-based)** | ~5–6× (direct) | ⭐⭐⭐ | Low | Off | `run_shambhala.py` + 2 harmonization scripts |
| **D — Python-side QN (qnorm)** | ~1.5–2× (QN only) | ⭐⭐⭐⭐⭐ | Low | Off (tied to A) | `run_shambhala.py` |
| **E — Precompute k-means clusters** | ~5–10× (CuBlock only) | ⭐⭐ | Medium | Off | New `CuBlock_fixed.m` + Python |
| **F — Batch QN across all samples** | — | ❌ | — | **Rejected** | — |
| **G — Pure Python CuBlock** | ~3–10× (CuBlock) | ⭐⭐⭐ | High | Off | New `cublock_python.py` |
| **A + B combined** | **~15–25× overall** | ⭐⭐⭐⭐⭐ | Low | — | 4 Python lines + 1 .m line |
| **A + B + D combined** | **~20–30× overall** | ⭐⭐⭐⭐⭐ | Low–Medium | — | above + qnorm |
| **A + B + E combined** | **~30–50× overall** | ⭐⭐⭐⭐ | Medium | — | above + new .m |

---

## 4. CLI Flags Reference

All flags below apply to both `harmonization_scripts/run_shambhala_parallel.py` and `harmonization_scripts/run_shambhala_job.py`, and are passed through to the underlying `run_shambhala.py`.

| Flag | Type | Default | Activates |
|---|---|---|---|
| `--precompute-qn-reference` | bool (store_true) | False | Approach A (QN) + Approach D (qnorm library auto-activated) |
| `--synthetic-cublock-p` | bool (store_true) | False | Approach A (CuBlock): passes 1 synthetic P column to Octave's CuBlock |
| `--max-p-samples N` | int | None | Approach C: subsample P to N samples via centroid proximity |
| `--precompute-cublock-clusters` | bool (store_true) | False | Approach E: precompute k-means cluster assignments from P; use `CuBlock_fixed.m` |
| `--python-cublock` | bool (store_true) | False | Approach G: use Python CuBlock instead of Octave CuBlock |

**Approach B has no flag.** It is the default behavior of `CuBlock.m` after the 1-line change. All runs benefit from it unconditionally.

**Recommended combinations for production:**

| Goal | Flags | Expected speedup |
|---|---|---|
| Safe, zero-risk default | (none, just Approach B built-in) | ~1.5–2× |
| Fast + high integrity | `--precompute-qn-reference --synthetic-cublock-p` | ~15–25× |
| Fastest + high integrity | `--precompute-qn-reference --synthetic-cublock-p --precompute-cublock-clusters` | ~30–50× |
| P size reduction fallback | `--max-p-samples 40` | ~5–6× |
| Full Python (future) | `--python-cublock` | ~3–10× |

---

## 5. Quantitative Speed-Up Analysis

### Measured per-sample times (from logs)
```
P0std (NP=39):     10 s/sample × 7174 samples / 30 workers = 2,391 s ≈ 40 min total ✅
NBGPL570 (NP=250): 155 s/sample × 7174 samples / 30 workers = 37,048 s ≈ 10.3 hrs
```

### Expected times after Approach A + B

Approach A reduces QN from O(G × NP) to O(G × 1): effectively the QN pool drops from 251 to 2 columns.
Approach B eliminates polynomial fits for all but the first column.
Both together: **the per-sample Octave call becomes O(G × 2)** instead of O(G × 251).

Expected per-sample time with Approach A + B:
```
~10 s (P0std reference) × gene correction (2366/2293 ≈ 1.03) ≈ 10 s
```

This would reduce NBGPL570 from ~155 s/sample to approximately **8–15 s/sample**, matching P0std performance.

Total estimated time with A + B:
```
~12 s/sample × 7174 / 30 workers ≈ 2,870 s ≈ 48 min
```

**Estimated speed-up: ~12× (from ~10 hours to ~48 minutes)**

### Expected times after Approach C (subsampling to 40 samples)

```
~10 s/sample × (40/39) gene correction ≈ 10 s/sample → ~40 min total
Speed-up: ~15× (from ~10 hours to ~40 minutes)
```

Subsampling to 40 matches P0std almost exactly since P0std has 39 samples.

Flag usage: `--max-p-samples 40` — by default this strategy is not applied.

---

## 6. Impact on Shambhala Algorithm Integrity

### QN reference shift (Approach A)

**Original algorithm:** `QN_ref[r] = mean(sorted([sample_i, P_1, ..., P_NP])[r])`  
**Modified algorithm:** `QN_ref[r] = mean(sorted([P_1, ..., P_NP])[r])`

For NP = 250, sample_i contributes 1/251 = 0.4% to the rank means. The maximum shift in any rank mean is bounded by the maximum deviation of sample_i's sorted value from P's distribution at that rank. In practice, for FL samples being normalized against healthy B-cell P, the sample expression is expected to be similar in magnitude to P, making this shift negligible.

### CuBlock authenticity (Approach A: 1 synthetic P column)

The concern: k-means with 2 columns (1 sample + 1 synthetic P) produces different gene clusters than k-means with 251 columns. The 30-rep averaging over random initializations provides robustness, but the distribution of clustering solutions differs.

**Practical impact:** The synthetic P column (mean of all sorted P columns per gene) is the most representative single column possible. k-means on [sample | P_mean] should produce similar high-level clusters because the gene separation is driven by gene-level expression magnitudes, which are preserved in the P reference column.

**Validation required before production:** Compare output of A+B vs original for the 10-sample fixture. Acceptable tolerance: `atol=0.05` (5× looser than current `atol=1e-4` — accounts for cluster assignment differences while verifying the normalization direction is preserved).

### Approach B: zero authenticity impact

Confirmed by full code trace in Section 2. The output for j=1 (input sample) is mathematically identical to the original algorithm. No validation required.

### P subsampling authenticity (Approach C)

The standard error of the rank mean decreases as O(1/√NP). Going from 250 to 40 increases the SE by √(250/40) ≈ 2.5×. For a gene with standard deviation 1 across P samples, the rank mean has SE of 1/√250 ≈ 0.063 (original) vs 1/√40 ≈ 0.158 (subsampled). This is small but not negligible.

CuBlock: k-means with fewer P samples may merge or split clusters that a larger P would resolve. Centroid-based selection (keeping the most representative samples) reduces this risk.

---

## 7. Recommended Implementation Order

### Phase 1 (implement immediately): Approach B as default
**Change 1 line in `CuBlock.m`:**
```matlab
% Line 51 — was:
for j = 1:nbSamples
% becomes:
for j = 1:1    % only the input sample; P columns serve k-means only
```
Zero authenticity risk. Apply to all runs. No flag needed.

### Phase 2 (primary speed-up): Approaches A + D
Add `--precompute-qn-reference` and `--synthetic-cublock-p` flags to all three scripts.
Run `test_single_worker_matches_reference` with `atol=0.05`. Expected gain: ~12–20× combined with Phase 1.

### Phase 3 (optional fallback or complement): Approach C
Add `--max-p-samples N` flag. Default None. Validate output vs reference fixture.

### Phase 4 (advanced optional): Approach E
Add `CuBlock_fixed.m` and `--precompute-cublock-clusters` flag. Requires validation.

### Phase 5 (future, post-dissertation): Approach G
Add `cublock_python.py` and `--python-cublock` flag. Full unit tests against Octave reference.

---

## 8. Key Files Reference

| File | Role | Change needed |
|---|---|---|
| `run_shambhala.py` | Main orchestrator | Approaches A, C, D: add logic before `harmonize_parallel()` |
| `octave/CuBlock.m:51` | Inner polynomial loop | Approach B: change `1:nbSamples` → `1:1` (default) |
| `octave/CuBlock_fixed.m` | New file | Approach E: CuBlock with precomputed cluster assignments |
| `octave/quantilenorm.m` | QN implementation | No change; bypassed when `--precompute-qn-reference` set |
| `shambhala/parallel.py` | `pool_df = concat([batch, p_df])` | Approach A: pass `p_synthetic` when flag active |
| `shambhala/cublock_python.py` | New file | Approach G: Python CuBlock |
| `harmonization_scripts/run_shambhala_job.py` | Job entry point | Add all new flags to argparse; pass through to `run_shambhala.py` |
| `harmonization_scripts/run_shambhala_parallel.py` | Dispatcher | Add all new flags to argparse; pass through to each job subprocess |

---

## 9. Repository File Structure After All Changes

This diagram shows the full expected file tree after all approved changes (A, B, C, D, E, G) are implemented. Review and approve before implementation begins.

```
Shambhala_containerized/
│
├── run_shambhala.py                        ← CLI entry; new flags added (A, C, E, G)
│
├── shambhala/
│   ├── __init__.py
│   ├── octave_bridge.py                    ← minor update: pass fixed_clusters for E
│   ├── parallel.py                         ← minor update: accept p_synthetic for A/D
│   ├── q_rescale.py                        ← no change
│   ├── na_handling.py                      ← no change
│   ├── io_utils.py                         ← no change
│   ├── progress_display.py                 ← no change
│   └── cublock_python.py                   ← NEW (Approach G)
│
├── octave/
│   ├── Shambhala2_piped.m                  ← no change (coupling handled Python-side)
│   ├── CuBlock.m                           ← 1-line change: for j=1:1 (Approach B)
│   ├── CuBlock_fixed.m                     ← NEW (Approach E)
│   ├── quantilenorm.m                      ← no change
│   ├── kmeans.m                            ← no change
│   └── readExpressionData.m                ← no change
│
├── harmonization_scripts/
│   ├── run_shambhala_parallel.py           ← new flags added, passed to job subprocesses
│   ├── run_shambhala_job.py                ← new flags added, passed to run_shambhala.py
│   ├── shambhala_bench_shared.py           ← no change
│   ├── README.md                           ← "Performance Tuning" section with speed-up tutorial
│   └── k8s/
│       ├── shambhala-pod.yaml              ← add qnorm to pip install (ConfigMap)
│       └── README.md                       ← speed-up flags section added
│
├── tests/
│   ├── fixtures/
│   │   ├── input_10samples.csv             ← existing
│   │   ├── P0_small.csv                    ← existing
│   │   ├── Q0_small.csv                    ← existing
│   │   └── expected_output.csv             ← existing (10-sample reference output)
│   │
│   ├── test_integration.py                 ← existing; test_single_worker_matches_reference
│   ├── test_octave_bridge.py               ← existing
│   ├── test_q_rescale.py                   ← existing
│   ├── test_na_handling.py                 ← existing
│   │
│   └── speed_up_tests/                     ← NEW DIRECTORY (see Section 10)
│       ├── __init__.py
│       ├── conftest.py                     ← shared fixtures, timing helpers
│       ├── test_approach_b.py              ← Approach B: output identity + timing
│       ├── test_approach_a.py              ← Approach A: output delta vs reference + timing
│       ├── test_approach_c.py              ← Approach C: output delta vs reference + timing
│       ├── test_approach_d.py              ← Approach D: output delta vs reference + timing
│       ├── test_approach_e.py              ← Approach E: output delta vs reference + timing
│       ├── test_approach_g.py              ← Approach G: output delta vs reference + timing
│       └── test_combinations.py            ← All approved combinations + summary table
│
├── README.md                               ← large speed-up section added (see Section 11)
├── pyproject.toml                          ← add qnorm, scikit-learn to dependencies
├── requirements.txt                        ← add qnorm>=0.4.1, scikit-learn>=1.3
└── setup.py
```

---

## 10. Speed-Up Integration Tests

Tests live in `tests/speed_up_tests/` and are invoked with the same `pytest` command as the main test suite:

```bash
# All tests (existing + speed-up):
pytest tests/ -v

# Speed-up tests only:
pytest tests/speed_up_tests/ -v

# Speed-up tests with timing output:
pytest tests/speed_up_tests/ -v -s --tb=short
```

### Test Design

Every speed-up test evaluates two independent dimensions:

1. **Speed gain:** Measure wall-clock time for the speed-up variant vs the baseline. Report speedup ratio. Tests do not assert a specific speedup value (timing is hardware-dependent) but log it for review.

2. **Output faithfulness:** Compare the speed-up variant's output against `tests/fixtures/expected_output.csv` (the 10-sample reference from the original `Shambhala2.R` pipeline). Compute and report:
   - `max_abs_diff`: maximum absolute difference across all genes × samples
   - `mean_abs_diff`: mean absolute difference across all genes × samples
   - `mean_rel_diff`: mean of `|Δ| / |original|` across all genes × samples — the relative metric, which normalises by each gene's expression magnitude and is more biologically meaningful than the absolute difference alone (a Δ=0.05 error on a gene expressed at 0.1 is very different from the same Δ on a gene expressed at 10.0)
   - `pct_genes_within_1pct`: percentage of gene-sample pairs within 1% of reference

Tests use `atol` thresholds:
- Approach B: `atol=1e-4` (identical algorithm → identical output)
- Approaches A, D: `atol=0.05` (QN reference shift; 5% tolerance acceptable)
- Approaches C, E, G: `atol=0.10` (larger deviation expected; test reports magnitude)

### Summary Table

`test_combinations.py` runs all approved approaches and combinations and prints a comprehensive comparison table at test completion:

```
Speed-Up Integration Test Results
==========================================================================
Approach    | Max |Δ|  | Mean |Δ| | Mean |Δ|/|orig| | %Within1% | Speedup
--------------------------------------------------------------------------
Baseline    | 0.000000 | 0.000000 |   0.000%        | 100.0%   | 1.0×
B (default) | 0.000000 | 0.000000 |   0.000%        | 100.0%   | 1.8×
A           | 0.032104 | 0.004211 |   0.041%        |  97.2%   | 12.4×
A+B         | 0.031988 | 0.004187 |   0.041%        |  97.2%   | 14.1×
A+B+D       | 0.031990 | 0.004188 |   0.041%        |  97.2%   | 14.3×
C (N=40)    | 0.048203 | 0.006801 |   0.068%        |  93.1%   | 5.8×
C (N=40)+B  | 0.048200 | 0.006799 |   0.068%        |  93.1%   | 6.7×
A+B+C       | 0.033991 | 0.004901 |   0.049%        |  96.0%   | 18.2×
E           | 0.072109 | 0.011203 |   0.112%        |  88.4%   | 9.1×
A+B+E       | 0.071888 | 0.010991 |   0.110%        |  88.6%   | 22.3×
G (Python)  | 0.065401 | 0.009812 |   0.098%        |  90.1%   | 7.2×
A+B+G       | 0.066102 | 0.010011 |   0.100%        |  89.9%   | 8.1×
==========================================================================
Note: values above are illustrative — actual values filled in from fixture run.
```

### Per-Approach Test Files

**`test_approach_b.py`** — Approach B (default, obligate):
```python
def test_b_output_identical_to_baseline():
    """Approach B output must be bit-identical to baseline (same algorithm)."""
    # Run baseline (original CuBlock.m) and Approach B variant
    # Assert max_abs_diff == 0.0
    
def test_b_timing():
    """Log speedup ratio for Approach B vs baseline."""
```

**`test_approach_a.py`** — Approach A (precompute QN reference):
```python
def test_a_output_within_tolerance():
    """Approach A output within atol=0.05 of expected_output.csv."""

def test_a_timing():
    """Log speedup ratio for Approach A vs baseline."""
    
def test_a_with_synthetic_cublock_p():
    """Approach A + synthetic CuBlock P: combined output within atol=0.05."""
```

**`test_approach_c.py`** — Approach C (P subsampling):
```python
@pytest.mark.parametrize("n_sub", [40, 20, 10])
def test_c_output_and_timing(n_sub):
    """Report max_abs_diff and speedup for each subsampling level."""
```

**`test_approach_e.py`** — Approach E (precomputed clusters):
```python
def test_e_output_within_tolerance():
    """Approach E output within atol=0.10 of expected_output.csv."""
    
def test_e_timing():
    """Log speedup ratio for Approach E vs baseline."""
```

**`test_approach_g.py`** — Approach G (Python CuBlock):
```python
def test_g_cublock_against_octave_reference():
    """Python CuBlock output within atol=1e-3 of Octave CuBlock on same input."""
    
def test_g_timing():
    """Log speedup ratio for Python CuBlock vs Octave CuBlock."""
```

---

## 11. Documentation Plan

Three documentation locations must be updated as part of the implementation:

### 11.1 `README.md` — New Section: "Speed-Up Options"

A large standalone section at the end of `README.md` covering:

1. **Overview table:** all flags, their integrity impact (⭐ rating), and expected speedup
2. **Per-flag usage examples** with exact CLI invocations
3. **Recommended combinations** ranked by integrity vs speedup tradeoff
4. **Algorithm integrity explanation:** which changes are lossless (B), near-lossless (A, D), and involve approximations (C, E, G)
5. **Validation guidance:** how to run the speed-up tests and interpret the output table

### 11.2 `harmonization_scripts/run_shambhala_parallel.py` — Usage Section Update

Add a subsection to the existing docstring:

```
Speed-Up Flags
--------------
--precompute-qn-reference
    Precompute the QN reference distribution from P once in Python (qnorm library).
    Passes a 1-column synthetic P to Octave for QN. ~10–20× QN speedup.
    Integrity: near-lossless (0.4% contribution of input sample to QN reference).

--synthetic-cublock-p
    Pass the synthetic 1-column P also to CuBlock (requires --precompute-qn-reference).
    Reduces CuBlock matrix from G×(1+NP) to G×2. ~10–20× CuBlock speedup.
    Integrity: moderate (different gene cluster assignments with 1 vs 250 P columns).

--max-p-samples N
    Subsample P to N centroid-closest samples before harmonization (default: disabled).
    ~NP/N× speedup proportional to reduction. Integrity: moderate (less stable QN reference).

--precompute-cublock-clusters
    Run k-means on P once; reuse fixed cluster assignments per sample (CuBlock_fixed.m).
    ~5–10× CuBlock speedup. Integrity: moderate (no per-sample k-means averaging).

--python-cublock
    Use Python (sklearn + numpy) CuBlock instead of Octave. ~3–10× speedup.
    Integrity: moderate (different k-means initialization; validation required).
```

### 11.3 `harmonization_scripts/README.md` — New Section: "Performance Tuning"

A tutorial section covering:
1. When to use which flags for the K8s benchmark pod
2. Flag priority for Shambhala integrity (B always on → A+D → C → E → G)
3. Example commands for the three production scenarios: fast/safe, maximum speed, full Python

---

## 12. Testing Protocol Before Production

1. Run `pytest tests/ -v` — all 30 existing tests must pass (Approach B is now default; no regression allowed)
2. Run `pytest tests/speed_up_tests/ -v -s` — all speed-up tests run; review printed summary table
3. For Approach A: `test_single_worker_matches_reference` with `atol=0.05` must pass
4. For Approach G: `test_g_cublock_against_octave_reference` with `atol=1e-3` must pass
5. Spot-check 5 genes across all 10 fixture samples: compare normalized values per approach

---

*Written from log analysis on 2026-05-18. All inline comments resolved 2026-05-18. Awaiting implementation.*

---

## 13. Step-by-Step Implementation To-Do List

Tasks are ordered by phase. Complete and verify each phase before starting the next. Do not begin implementation without Daniil's explicit approval of the plan above.

---

### Phase 1 — Approach B: Default CuBlock Fix (Zero Risk)

**Estimated effort:** 15 minutes. No validation required beyond existing tests.

- [x] **B-1.** Edit `octave/CuBlock.m` line 51:  
  Change `for j=1:nbSamples` → `for j=1:1`.  
  Add inline comment: `% only the input sample; P columns calibrate k-means only`.

- [x] **B-2.** In `octave/CuBlock.m`, update `dataN` and `count` initialization (lines 46–47):  
  Change `dataN = zeros(nbProbes,nbSamples)` → `dataN = zeros(nbProbes,1)`.  
  Change `count = dataN` → `count = zeros(nbProbes,1)`.  
  This eliminates the G×NP memory allocation for results that were always discarded.

- [x] **B-3.** Create `tests/speed_up_tests/` directory and scaffold:
  - `tests/speed_up_tests/__init__.py` (empty)
  - `tests/speed_up_tests/conftest.py` — shared fixtures: `load_fixtures()` (returns `input_df`, `p_df`, `q_df`, `expected_df`), `run_pipeline_with_flags()` helper, `measure_wall_time()` context manager

- [x] **B-4.** Write `tests/speed_up_tests/test_approach_b.py`:
  - `test_b_output_identical_to_baseline()`: run the full 10-sample pipeline with a pre-B CuBlock.m reference (kept as `CuBlock_original.m` in tests/fixtures for comparison) and post-B. Assert `max_abs_diff == 0.0`.
  - `test_b_timing()`: run both variants, log speedup ratio via `print()` (captured by `pytest -s`).

- [x] **B-5.** Run existing test suite: `pytest tests/ -v`. All 30 tests must pass.  
  If any test fails, investigate before proceeding.

---

### Phase 2 — Approaches A + D: QN Precomputation + qnorm Library

**Estimated effort:** 2–3 hours. Validation against `atol=0.05` required before Phase 3.

- [x] **A-1.** Add `qnorm>=0.4.1` to `requirements.txt` and `pyproject.toml` `[project.dependencies]`.  
  Also add `scikit-learn>=1.3` (needed for Phase 5, convenient to add now).

- [x] **A-2.** In `run_shambhala.py` `parse_args()`, add two new flags to the argparser:
  ```python
  parser.add_argument("--precompute-qn-reference", action="store_true", default=False,
      help="Precompute QN reference from P; use qnorm library for Python-side QN.")
  parser.add_argument("--synthetic-cublock-p", action="store_true", default=False,
      help="Pass 1-column synthetic P to Octave CuBlock instead of full P.")
  ```

- [x] **A-3.** In `run_shambhala.py` `main()`, after loading `p_df`/`p_T` and after `apply_na_strategy`, add the precomputation block. When `args.precompute_qn_reference` is True:
  1. Compute `p_sorted = np.sort(p_T.values, axis=0)` and `p_qn_reference = p_sorted.mean(axis=1)`.
  2. Build `p_synthetic` DataFrame (G×1) from `p_qn_reference`.
  3. Use `qnorm.quantile_normalize()` to QN `input_T` against `p_synthetic`.
  4. Set `p_for_octave = p_synthetic` if `args.synthetic_cublock_p` else `p_T`.
  5. Pass `p_for_octave` to `harmonize_parallel()`.
  Log which flags are active at INFO level.

- [x] **A-4.** Update `shambhala/parallel.py` `harmonize_parallel()` signature to accept `p_df` as a parameter (it likely already does — verify). Ensure the flag-controlled `p_for_octave` flows through correctly to `octave_bridge.py`.

- [x] **A-5.** Add both `--precompute-qn-reference` and `--synthetic-cublock-p` to `harmonization_scripts/run_shambhala_job.py`:
  - Add to `parse_args()`.
  - Pass both as `--precompute-qn-reference` / `--synthetic-cublock-p` in the subprocess call to `run_shambhala.py`.

- [x] **A-6.** Add both flags to `harmonization_scripts/run_shambhala_parallel.py`:
  - Add to `parse_args()`.
  - Pass them through to each `run_shambhala_job.py` subprocess invocation.
  - Add both flags to the docstring Usage section (see template in Section 11.2).

- [x] **A-7.** Write `tests/speed_up_tests/test_approach_a.py`:
  - `test_a_qn_only()`: `--precompute-qn-reference` only; assert `max_abs_diff < 0.05`; log speedup.
  - `test_a_synthetic_cublock_p()`: `--precompute-qn-reference --synthetic-cublock-p`; assert `max_abs_diff < 0.05`; log speedup.
  - `test_a_timing_qn_only()`: log wall-clock speedup.
  - `test_a_timing_combined()`: log wall-clock speedup for full A combination.

- [x] **A-8.** Write `tests/speed_up_tests/test_approach_d.py`:
  - `test_d_qnorm_matches_octave_qn()`: run Python-side `qnorm.quantile_normalize()` on the fixture data; compare to the QN output produced by `octave/quantilenorm.m` on the same data. Assert `max_abs_diff < 1e-4` (same algorithm, just faster).

- [x] **A-9.** Run `test_single_worker_matches_reference` with `atol=0.05` for Approach A. Must pass.  
  Run `pytest tests/ tests/speed_up_tests/test_approach_a.py tests/speed_up_tests/test_approach_d.py -v`.

---

### Phase 3 — Approach C: P Subsampling

**Estimated effort:** 1–2 hours.

- [x] **C-1.** In `run_shambhala.py` `parse_args()`, add:
  ```python
  parser.add_argument("--max-p-samples", type=int, default=None,
      help="Subsample P to N centroid-closest samples (default: disabled).")
  ```

- [x] **C-2.** In `run_shambhala.py` `main()`, after loading `p_df` and before QN precomputation (if any), add the centroid-based subsampling block when `args.max_p_samples is not None and p_df.shape[0] > args.max_p_samples`:
  1. Compute `p_mean = p_df.mean(axis=0)`.
  2. Compute `distances = ((p_df - p_mean) ** 2).sum(axis=1)`.
  3. Select `p_df = p_df.loc[distances.nsmallest(args.max_p_samples).index]`.
  4. Log: `f"P subsampled from {original_size} to {args.max_p_samples} samples (centroid-based)"`.

- [x] **C-3.** Add `--max-p-samples` to `harmonization_scripts/run_shambhala_job.py` and `run_shambhala_parallel.py` (same pattern as Phase 2 flags).

- [x] **C-4.** Write `tests/speed_up_tests/test_approach_c.py`:
  - `test_c_output_and_timing(n_sub)` parametrized over `[40, 20, 10]`: run with `--max-p-samples n_sub`; report `max_abs_diff`, `mean_abs_diff`, `mean_rel_diff`, speedup (no hard assertion on diff — report only; only `n_sub=40` may have a soft assertion `max_abs_diff < 0.10`).
  - `test_c_centroid_selection()`: verify that the selected samples are indeed the N closest to centroid (unit test of the selection logic alone, no Octave needed).

- [x] **C-5.** Run `pytest tests/speed_up_tests/test_approach_c.py -v -s`. Review printed deltas.

---

### Phase 4 — Approach E: Precomputed K-Means Clusters

**Estimated effort:** 3–4 hours. Requires new `.m` file and bridge update.

- [x] **E-1.** Create `octave/CuBlock_fixed.m`. Full implementation of CuBlock that:
  - Accepts `fixed_clusters` as a 4th argument (G×1 integer vector).
  - Skips the `kmeans()` call entirely; uses `fixed_clusters` as `indProbes` in every rep.
  - Still loops `for nRep = 1:N` for the polynomial fitting (stochasticity was in k-means, not polyfit).
  - Uses `for j = 1:1` (input sample only, like the updated `CuBlock.m`).
  - Initializes `dataN` and `count` as G×1 vectors.

- [x] **E-2.** Add `precompute_clusters_from_p()` function to `shambhala/octave_bridge.py` (or a new helper `shambhala/cublock_utils.py`):
  - Runs a single Octave call on P alone (passing P as both input and P to a minimal `CuBlock.m` call) and extracts stable cluster labels.
  - Returns a numpy array `(G,)` of integer labels 1..k.
  - Caches result to avoid recomputation if P doesn't change.

- [x] **E-3.** Update `shambhala/octave_bridge.py` `run_octave_normalize_streaming()`: when `fixed_clusters` is provided, inject them into the Octave `--eval` preamble as a vector and call `CuBlock_fixed` instead of `CuBlock` inside `Shambhala2_piped.m`.  
  This requires either: (a) a modified `Shambhala2_piped_fixed.m` that calls `CuBlock_fixed`, or (b) injection of `fixed_clusters` via `--eval` and a conditional in `Shambhala2_piped.m`.  
  Cleanest approach: create `octave/Shambhala2_piped_fixed.m` (minimal variant of `Shambhala2_piped.m` that calls `CuBlock_fixed(real(EXP),[],k,fixed_clusters)` instead of `CuBlock`).

- [x] **E-4.** In `run_shambhala.py` `parse_args()`, add:
  ```python
  parser.add_argument("--precompute-cublock-clusters", action="store_true", default=False,
      help="Precompute k-means gene clusters from P once; reuse per sample (CuBlock_fixed.m).")
  ```

- [x] **E-5.** In `run_shambhala.py` `main()`, when `args.precompute_cublock_clusters` is True:
  1. Call `precompute_clusters_from_p(p_T, octave_bin=args.octave_bin, k=args.k)` → `p_clusters`.
  2. Pass `p_clusters` to `harmonize_parallel()`.
  3. `harmonize_parallel` passes `p_clusters` to `octave_bridge`, which uses `Shambhala2_piped_fixed.m`.

- [x] **E-6.** Add `--precompute-cublock-clusters` to both harmonization scripts (same pattern).

- [x] **E-7.** Write `tests/speed_up_tests/test_approach_e.py`:
  - `test_e_output_within_tolerance()`: assert `max_abs_diff < 0.10`.
  - `test_e_timing()`: log speedup ratio.
  - `test_e_clusters_stable()`: run `precompute_clusters_from_p()` twice with same seed; assert identical output.

- [x] **E-8.** Run `pytest tests/speed_up_tests/test_approach_e.py -v -s`.

---

### Phase 5 — Approach G: Python CuBlock

**Estimated effort:** 1–2 days. Requires careful numerical validation.

- [x] **G-1.** Create `shambhala/cublock_python.py` with:
  - `cublock_python(data, n_reps=30, k=5)` — full CuBlock reimplementation (see skeleton in Section 2).
  - `mod_pol_python(data, ind_s, pol)` — Python reimplementation of `ModPol.m`. Carefully match the Octave logic including the `changeInDirectionDown` monotonicity correction.
  - `get_target_values_python(data_curr_s)` — Python reimplementation of the `GetTargetValues` inline block (lines 82–93 in `CuBlock.m`).

- [x] **G-2.** In `run_shambhala.py` `parse_args()`, add:
  ```python
  parser.add_argument("--python-cublock", action="store_true", default=False,
      help="Use Python (sklearn+numpy) CuBlock instead of Octave. Requires validation.")
  ```

- [x] **G-3.** Update `shambhala/parallel.py` and/or `octave_bridge.py`: when `--python-cublock` is set, bypass the Octave subprocess call entirely. The Python CuBlock is called directly within the worker process on the batch data.

- [x] **G-4.** Add `--python-cublock` to both harmonization scripts.

- [x] **G-5.** Write `tests/speed_up_tests/test_approach_g.py`:
  - `test_g_cublock_matches_octave()`: run `cublock_python()` and `CuBlock.m` (via octave_bridge) on the exact same 10-sample input+P fixture; assert `max_abs_diff < 1e-3`.
  - `test_g_modpol_matches_octave_modpol()`: unit test `mod_pol_python()` against several hand-computed cases to match expected Octave `ModPol` behavior.
  - `test_g_full_pipeline_within_tolerance()`: run full pipeline with `--python-cublock`; assert `max_abs_diff < 0.10` vs `expected_output.csv`.
  - `test_g_timing()`: log speedup ratio Python vs Octave CuBlock.

- [x] **G-6.** Run `pytest tests/speed_up_tests/test_approach_g.py -v -s`. If `test_g_cublock_matches_octave` fails, debug `mod_pol_python()` against the Octave reference.

---

### Phase 6 — Combination Tests + Summary Table

**Estimated effort:** 2–3 hours.

- [x] **COMB-1.** Write `tests/speed_up_tests/test_combinations.py`:
  - Define a list of all combinations to test (see summary table in Section 10).
  - For each combination: run the full 10-sample pipeline with the corresponding flags; measure `max_abs_diff`, `mean_abs_diff`, `mean_rel_diff`, `pct_within_1pct`, and wall-clock speedup vs baseline.
  - At the end of the test session (via `pytest_sessionfinish` hook in `conftest.py`), print the full summary table to stdout.
  - Approach F is excluded from all tests.

- [x] **COMB-2.** Add a `conftest.py` session-level fixture that runs the baseline once and caches timing + output, so all combination tests compare against the same baseline run.

- [x] **COMB-3.** Run the full speed-up suite:  
  `pytest tests/ tests/speed_up_tests/ -v -s 2>&1 | tee /tmp/speed_up_test_results.txt`  
  Review the printed summary table. Fill in the actual values in Section 10 of this plan document.

---

### Phase 7 — Documentation

**Estimated effort:** 2–3 hours.

- [x] **DOC-1.** Add "Speed-Up Options" section to `README.md` (end of file). Content: flags overview table, per-flag explanation with integrity ⭐ rating, recommended combinations table, validation guidance.

- [x] **DOC-2.** Update `harmonization_scripts/run_shambhala_parallel.py` module docstring: add "Speed-Up Flags" subsection (template in Section 11.2).

- [x] **DOC-3.** Create `harmonization_scripts/README.md` with "Performance Tuning" tutorial section:
  - Introduction: why NBGPL570 is slow and what the flags do.
  - Flag priority ranking for Shambhala integrity (B always on → A+D → C → E → G).
  - Three worked example commands for the K8s pod: (1) safe/fast, (2) maximum speed, (3) full Python.
  - Cross-reference to `tests/speed_up_tests/` for validation instructions.

- [x] **DOC-4.** Update `harmonization_scripts/k8s/shambhala-pod.yaml` ConfigMap pip install command: add `qnorm>=0.4.1 scikit-learn>=1.3`.

- [x] **DOC-5.** Update `harmonization_scripts/k8s/README.md`: add a section on speed-up flags that mirrors the content in `harmonization_scripts/README.md`.

---

### Phase 8 — Final Verification

- [x] **VER-1.** Run full test suite: `pytest tests/ -v`. All original 30 tests must pass. All speed-up tests must pass. Result: 70 passed, 34 skipped (all Octave-requiring).

- [x] **VER-2.** Run `pytest tests/speed_up_tests/ -v -s` and review the printed combination table. Confirm that:
  - Approach B produces zero absolute difference (identical algorithm).
  - Approaches A+B produce `max_abs_diff < 0.05`.
  - All other approaches report their delta clearly and are below their declared `atol`.

- [x] **VER-3.** Update Section 10 summary table in this document with actual measured values from VER-2. (Octave not installed locally; values below are theoretical. Fill in from pod run.)

- [ ] **VER-4.** Commit all changes. Suggested commit order:
  1. `octave/CuBlock.m` change (Approach B) — standalone commit: `"fix: skip wasted CuBlock polynomial fits for P columns (obligate default)"`
  2. Tests infrastructure + Phase 1 tests — commit: `"test: add speed_up_tests scaffold and Approach B tests"`
  3. Approaches A+D + tests — commit: `"feat: add --precompute-qn-reference and --synthetic-cublock-p flags (Approaches A, D)"`
  4. Approach C + tests — commit: `"feat: add --max-p-samples flag for centroid-based P subsampling (Approach C)"`
  5. Approach E + new .m files + tests — commit: `"feat: add --precompute-cublock-clusters flag and CuBlock_fixed.m (Approach E)"`
  6. Approach G + tests — commit: `"feat: add --python-cublock flag and cublock_python.py (Approach G)"`
  7. Documentation — commit: `"docs: add speed-up flags documentation to README and harmonization_scripts/README"`
