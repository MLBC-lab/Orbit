from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    n = len(p)
    adjusted = np.empty(n, dtype=float)
    cumulative = 1.0
    for i in range(n - 1, -1, -1):
        cumulative = min(cumulative, ranked[i] * n / (i + 1))
        adjusted[order[i]] = cumulative
    return np.clip(adjusted, 0.0, 1.0)


def expression_to_signature(
    expression: pd.DataFrame,
    samples: pd.DataFrame,
    *,
    sample_col: str = "sample",
    group_col: str = "group",
    case_label: str = "case",
    control_label: str = "control",
) -> pd.DataFrame:
    required = {sample_col, group_col}
    if not required.issubset(samples.columns):
        raise ValueError(f"samples table must contain {sorted(required)}")
    sample_order = samples[sample_col].astype(str).tolist()
    missing = [s for s in sample_order if s not in expression.columns]
    if missing:
        raise ValueError(f"expression table is missing sample columns: {missing[:5]}")
    expr = expression.set_index("gene")[sample_order]
    groups = samples.set_index(sample_col).loc[sample_order, group_col]
    case_cols = groups[groups == case_label].index.tolist()
    control_cols = groups[groups == control_label].index.tolist()
    if len(case_cols) < 2 or len(control_cols) < 2:
        raise ValueError("at least two case and two control samples are required")
    case = expr[case_cols].astype(float)
    control = expr[control_cols].astype(float)
    t_stat, p_value = stats.ttest_ind(case.T, control.T, equal_var=False, nan_policy="omit")
    effect = case.mean(axis=1).to_numpy() - control.mean(axis=1).to_numpy()
    out = pd.DataFrame(
        {
            "gene": expr.index.astype(str),
            "effect": effect,
            "statistic": np.nan_to_num(t_stat),
            "p_value": np.nan_to_num(p_value, nan=1.0),
        }
    )
    out["q_value"] = benjamini_hochberg(out["p_value"].to_numpy())
    return out


def topn_sets(signature: pd.DataFrame, n: int) -> tuple[set[str], set[str]]:
    ranked = signature.sort_values("effect", ascending=False)
    up = set(ranked.head(n)["gene"].astype(str))
    down = set(ranked.tail(n)["gene"].astype(str))
    return up, down


def signature_vector(signature: pd.DataFrame, genes: list[str]) -> np.ndarray:
    values = signature.groupby("gene")["effect"].mean()
    return values.reindex(genes).fillna(0.0).to_numpy(dtype=float)


def stouffer_consensus(signatures: pd.DataFrame, weight_col: str | None = None) -> pd.DataFrame:
    records = []
    for (disease, method, gene), group in signatures.groupby(["disease", "method", "gene"]):
        weights = np.ones(len(group)) if weight_col is None else group[weight_col].to_numpy(float)
        z = group["effect"].to_numpy(float)
        denom = np.sqrt(np.sum(weights**2))
        combined = float(np.sum(weights * z) / denom) if denom else 0.0
        p_value = float(2.0 * stats.norm.sf(abs(combined)))
        records.append(
            {
                "disease": disease,
                "cohort": "consensus_stouffer",
                "method": method,
                "gene": gene,
                "effect": combined,
                "p_value": p_value,
            }
        )
    out = pd.DataFrame.from_records(records)
    out["q_value"] = out.groupby(["disease", "method"])["p_value"].transform(
        lambda x: benjamini_hochberg(x.to_numpy())
    )
    return out


def rank_aggregation_consensus(signatures: pd.DataFrame) -> pd.DataFrame:
    records = []
    for (disease, method), block in signatures.groupby(["disease", "method"]):
        pivot = block.assign(
            rank=block.groupby("cohort")["effect"].rank(ascending=False, method="average")
        ).pivot_table(index="gene", columns="cohort", values="rank", aggfunc="mean")
        median_rank = pivot.median(axis=1)
        centered = -(median_rank - median_rank.median())
        for gene, effect in centered.items():
            records.append(
                {
                    "disease": disease,
                    "cohort": "consensus_rra",
                    "method": method,
                    "gene": gene,
                    "effect": float(effect),
                    "p_value": 1.0,
                    "q_value": 1.0,
                }
            )
    return pd.DataFrame.from_records(records)

