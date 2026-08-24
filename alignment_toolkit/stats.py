
import numpy as np
from scipy.stats import norm

from config import MIN_N_FOR_CI

def median_ci(values, ci):
    v = np.asarray(values, dtype=float)
    v = np.sort(v[~np.isnan(v)])
    n = v.size

    if n == 0:
        return np.nan, np.nan, np.nan
    median = float(np.median(v))
    if n < MIN_N_FOR_CI:
        return median, np.nan, np.nan

    z = norm.ppf(0.5 + ci / 200.0)          # z-score for the CI level
    half_width = z * np.sqrt(n) / 2.0        
    low_rank = max(int(np.floor(n / 2.0 - half_width)), 0)
    high_rank = min(int(np.ceil(n / 2.0 + half_width)) - 1, n - 1)

    return median, float(v[low_rank]), float(v[high_rank])


def quartiles(values):
    v = np.asarray(values, dtype=float)
    v = v[~np.isnan(v)]
    if v.size == 0:
        return np.nan, np.nan
    q1, q3 = np.percentile(v, [25, 75])
    return float(q1), float(q3)