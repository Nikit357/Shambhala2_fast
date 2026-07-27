"""
Pure Python reimplementation of CuBlock.m (Approach G).

Replaces the Octave subprocess for CuBlock normalization. The algorithm is
identical to CuBlock.m (as modified by Approach B: j=1:1 only). K-means uses
sklearn with random init and n_init=1 per rep, matching the Octave custom
kmeans.m pure-random initialization.
"""

import numpy as np
from sklearn.cluster import KMeans

_P_EXPONENTS = np.array([3, 5, 7, 9, 11, 13, 15, 17, 19, 21])
_TOL = 0.1
_MIN_CLUSTER_SIZE = 100


def mod_pol_python(data: np.ndarray, ind_s: np.ndarray, pol: np.ndarray) -> np.ndarray:
    """
    Evaluate cubic polynomial on sorted data and apply monotonicity correction.

    Direct translation of ModPol.m. Returns a vector aligned with the original
    (unsorted) `data` ordering; positions not covered by ind_s receive NaN.

    Parameters
    ----------
    data : np.ndarray shape (n,)
        Z-transformed block values in original (unsorted) order.
    ind_s : np.ndarray shape (n,) int
        Argsort of `data` (0-indexed).
    pol : np.ndarray shape (4,)
        Cubic polynomial coefficients [a3, a2, a1, a0].

    Returns
    -------
    np.ndarray shape (n,)
    """
    x_val = data[ind_s]
    data_n_s = pol[0] * x_val**3 + pol[1] * x_val**2 + pol[2] * x_val + pol[3]
    n = len(data)

    diffs = data_n_s[1:] - data_n_s[:-1]
    down_pos = np.where(diffs < 0)[0]

    if len(down_pos) > 0:
        # Octave 1-indexed positions (for slice arithmetic matching Octave exactly)
        ind_down1_oct = down_pos[0] + 1          # indDown1
        ind_down_l_oct = down_pos[-1] + 2         # find(...,'last') + 1

        up_pos = np.where(diffs > 0)[0]
        if len(up_pos) > 0:
            ind_up1_oct = up_pos[0] + 1            # indUp1
            ind_up_l_oct = up_pos[-1] + 2          # find(...,'last') + 1

            # 0-indexed element access (Octave 1-indexed - 1)
            idx_up1 = ind_up1_oct - 1
            idx_up_l = ind_up_l_oct - 1

            if data_n_s[idx_up1] == data_n_s[0] and data_n_s[idx_up_l] == data_n_s[n - 1]:
                big_m = data_n_s[ind_down1_oct - 1]   # M = dataNS(indDown1)
                small_m = data_n_s[ind_down_l_oct - 1] # m = dataNS(indDownL)
                lo = ind_down1_oct - 1                  # 0-indexed start of bump
                hi = ind_down_l_oct                     # 0-indexed exclusive end
                if (n - ind_down_l_oct + 1) <= ind_down1_oct:
                    data_n_s[lo:hi] = big_m
                    data_n_s[hi:] = big_m + data_n_s[hi:] - small_m
                else:
                    data_n_s[lo:hi] = small_m
                    data_n_s[:lo] = small_m + data_n_s[:lo] - big_m
            else:
                big_m = data_n_s[idx_up_l]
                small_m = data_n_s[idx_up1]
                data_n_s[idx_up_l:] = big_m
                data_n_s[: idx_up1 + 1] = small_m

    result = np.full(n, np.nan)
    result[ind_s] = data_n_s
    return result


def cublock_python(
    data: np.ndarray,
    n_reps: int = 30,
    k: int = 5,
    random_seed: int | None = None,
) -> np.ndarray:
    """
    Python reimplementation of CuBlock.m with Approach-B default (j=1 only).

    K-means is run on the full data matrix (all columns) to determine gene
    clusters, but polynomial fitting is applied only to column 0 (the input
    sample). Results for all 30 reps are averaged per gene.

    Parameters
    ----------
    data : np.ndarray shape (n_genes, n_samples)
        QN-normalized log2 expression matrix. First column = input sample.
        All columns used for k-means clustering. Must be NaN-free.
    n_reps : int
        CuBlock repetitions. Default 30.
    k : int
        Gene clusters for k-means. Default 5.
    random_seed : int or None
        Seed for reproducible k-means. Each rep gets a derived sub-seed.

    Returns
    -------
    np.ndarray shape (n_genes, 1)
        CuBlock-normalized first column. NaN for genes whose cluster never
        reached 100 members or whose Z-transform std was always zero.
    """
    data = np.asarray(data, dtype=float)
    n_genes = data.shape[0]
    data_n = np.zeros(n_genes)
    count = np.zeros(n_genes)

    rng = np.random.default_rng(random_seed)
    data_col = data[:, 0]

    for _ in range(n_reps):
        rep_seed = int(rng.integers(0, 2**31)) if random_seed is not None else None
        km = KMeans(
            n_clusters=k,
            n_init=1,
            max_iter=100,
            init="random",
            random_state=rep_seed,
        )
        labels = km.fit_predict(data)

        for cluster_id in range(k):
            mask = labels == cluster_id
            if mask.sum() <= _MIN_CLUSTER_SIZE:
                continue
            data_curr = data_col[mask].copy()
            valid = ~np.isnan(data_curr)
            n_valid = int(valid.sum())
            if n_valid < 2:
                continue
            curr_mean = float(data_curr[valid].mean())
            curr_std = float(data_curr[valid].std(ddof=1))
            if curr_std <= 0.0:
                continue

            data_curr = (data_curr - curr_mean) / curr_std
            ind_s = np.argsort(data_curr)
            data_curr_s = data_curr[ind_s]
            n_curr = len(data_curr)

            # GetTargetValues: find optimal polynomial exponent
            X = np.linspace(-1.0, 1.0, n_curr)[:, np.newaxis] ** _P_EXPONENTS
            ind_std_up = int(np.argmin(np.abs(data_curr_s - 1.0)))
            ind_std_down = int(np.argmin(np.abs(data_curr_s + 1.0)))
            lo = min(ind_std_down, ind_std_up)
            hi = max(ind_std_down, ind_std_up)
            X_subset = np.abs(X[lo : hi + 1, :])
            if X_subset.shape[0] == 0:
                X_subset = np.abs(X[[lo], :])
            S = X_subset.mean(axis=0)
            below = np.where(S < _TOL)[0]
            ind_p = int(below[0]) if len(below) > 0 else len(_P_EXPONENTS) - 1

            # Cubic polyfit: solve V @ pol = X[:, ind_p]
            V = np.column_stack(
                [data_curr_s**3, data_curr_s**2, data_curr_s, np.ones(n_curr)]
            )
            pol, _, _, _ = np.linalg.lstsq(V, X[:, ind_p], rcond=None)

            curr_data_n = mod_pol_python(data_curr, ind_s, pol)
            data_n[mask] += curr_data_n
            count[mask] += 1

    with np.errstate(invalid="ignore", divide="ignore"):
        result = np.where(count > 0, data_n / count, np.nan)
    return result[:, np.newaxis]
