from __future__ import annotations

from collections.abc import Callable

import numpy as np
import pandas as pd
from scipy import stats


def common_score_vectors(a: pd.DataFrame, b: pd.DataFrame, score_col: str = "score") -> tuple[np.ndarray, np.ndarray]:
    merged = a[["perturbagen", score_col]].merge(
        b[["perturbagen", score_col]],
        on="perturbagen",
        suffixes=("_a", "_b"),
    )
    if len(merged) < 3:
        return np.array([]), np.array([])
    return merged[f"{score_col}_a"].to_numpy(float), merged[f"{score_col}_b"].to_numpy(float)


def kendall_tau(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 3 or len(b) < 3:
        return np.nan
    return float(stats.kendalltau(a, b, nan_policy="omit").statistic)


def spearman_rho(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 3 or len(b) < 3:
        return np.nan
    return float(stats.spearmanr(a, b, nan_policy="omit").statistic)


def topk_set(frame: pd.DataFrame, k: int, score_col: str = "score") -> set[str]:
    return set(frame.sort_values(score_col, ascending=False).head(k)["perturbagen"].astype(str))


def jaccard_at_k(a: pd.DataFrame, b: pd.DataFrame, k: int, score_col: str = "score") -> float:
    left = topk_set(a, k, score_col)
    right = topk_set(b, k, score_col)
    union = left | right
    if not union:
        return np.nan
    return len(left & right) / len(union)


def overlap_at_k(a: pd.DataFrame, b: pd.DataFrame, k: int, score_col: str = "score") -> float:
    left = topk_set(a, k, score_col)
    right = topk_set(b, k, score_col)
    denom = min(len(left), len(right))
    if denom == 0:
        return np.nan
    return len(left & right) / denom


def rank_stability(a: pd.DataFrame, b: pd.DataFrame, top_k: tuple[int, ...], score_col: str = "score") -> dict[str, float]:
    x, y = common_score_vectors(a, b, score_col)
    out = {
        "n_common": int(len(x)),
        "kendall_tau": kendall_tau(x, y),
        "spearman_rho": spearman_rho(x, y),
    }
    for k in top_k:
        out[f"jaccard_at_{k}"] = jaccard_at_k(a, b, k, score_col)
        out[f"overlap_at_{k}"] = overlap_at_k(a, b, k, score_col)
    return out


def percentile_ci(values: np.ndarray, alpha: float = 0.05) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return np.nan, np.nan
    return (
        float(np.quantile(values, alpha / 2.0)),
        float(np.quantile(values, 1.0 - alpha / 2.0)),
    )


def bootstrap_metric(
    items: pd.DataFrame,
    metric: Callable[[pd.DataFrame], float],
    *,
    iterations: int,
    seed: int,
) -> tuple[float, float]:
    if len(items) == 0 or iterations <= 0:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(iterations):
        take = rng.integers(0, len(items), size=len(items))
        values.append(metric(items.iloc[take].reset_index(drop=True)))
    return percentile_ci(np.asarray(values))

