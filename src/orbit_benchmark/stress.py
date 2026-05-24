from __future__ import annotations

import numpy as np
import pandas as pd


def identify_stress_programs(
    loadings: pd.DataFrame,
    stress_gene_sets: pd.DataFrame,
    *,
    min_gene_hits: int = 5,
    top_fraction: float = 0.1,
) -> pd.DataFrame:
    stress_genes = set(stress_gene_sets["gene"].astype(str))
    records = []
    for program, block in loadings.groupby("program"):
        ranked = block.assign(abs_loading=block["loading"].abs()).sort_values("abs_loading", ascending=False)
        n_top = max(min_gene_hits, int(np.ceil(len(ranked) * top_fraction)))
        top = set(ranked.head(n_top)["gene"].astype(str))
        hits = len(top & stress_genes)
        expected = n_top * (len(stress_genes) / max(1, block["gene"].nunique()))
        enrichment = hits / max(expected, 1e-9)
        records.append(
            {
                "program": program,
                "stress_gene_hits": hits,
                "expected_hits": expected,
                "enrichment": enrichment,
                "is_stress_program": bool(hits >= min_gene_hits and enrichment >= 1.5),
            }
        )
    result = pd.DataFrame.from_records(records)
    if result["is_stress_program"].sum() == 0 and len(result):
        idx = result["enrichment"].idxmax()
        result.loc[idx, "is_stress_program"] = True
    return result


def stress_scores(
    activations: pd.DataFrame,
    stress_programs: pd.DataFrame,
    *,
    pooling: str = "sum",
) -> pd.Series:
    selected = stress_programs.loc[stress_programs["is_stress_program"], "program"].tolist()
    selected = [p for p in selected if p in activations.columns]
    if not selected:
        return pd.Series(0.0, index=activations.index, name="stress_score")
    values = activations[selected]
    if pooling == "max":
        score = values.abs().max(axis=1)
    else:
        score = values.sum(axis=1)
    return score.rename("stress_score")


def zscore_within_strata(frame: pd.DataFrame, group_cols: list[str], value_col: str) -> pd.Series:
    def z(block: pd.Series) -> pd.Series:
        std = block.std(ddof=0)
        if std == 0 or np.isnan(std):
            return pd.Series(0.0, index=block.index)
        return (block - block.mean()) / std

    return frame.groupby(group_cols, dropna=False)[value_col].transform(z)


def apply_stress_penalty(rankings: pd.DataFrame, stress_lambda: float) -> pd.DataFrame:
    out = rankings.copy()
    out["adjusted_score"] = out["score"] - stress_lambda * out["stress_z"]
    return out

