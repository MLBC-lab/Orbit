from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from .config import OrbitConfig
from .data import OrbitTables, pivot_perturbations
from .fair import write_manifest, write_ro_crate
from .metrics import rank_stability
from .objectives import (
    connectivity_reversal,
    learn_program_basis,
    program_reversal,
    vector_reversal,
)
from .signatures import signature_vector, topn_sets
from .stress import apply_stress_penalty, identify_stress_programs, stress_scores, zscore_within_strata


def quality_tier(tas: float, thresholds: dict[str, float]) -> str:
    ordered = sorted(thresholds.items(), key=lambda item: item[1], reverse=True)
    for name, threshold in ordered:
        if tas >= threshold:
            return name
    return "low"


def _disease_signature_groups(disease_signatures: pd.DataFrame):
    for keys, block in disease_signatures.groupby(["disease", "cohort", "method"], sort=True):
        yield keys, block.sort_values("gene")


def make_rankings(tables: OrbitTables, config: OrbitConfig) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    perturb_matrix, _ = pivot_perturbations(tables.perturbation_signatures)
    genes = sorted(set(tables.disease_signatures["gene"]).intersection(perturb_matrix.columns.astype(str)))
    if not genes:
        raise ValueError("no shared genes between disease and perturbation signatures")
    perturb_matrix = perturb_matrix.reindex(columns=genes).fillna(0.0)
    basis = learn_program_basis(perturb_matrix, config.program_k, config.random_seed)
    perturb_programs = basis.transform_matrix(perturb_matrix)
    loadings = basis.loadings()
    stress_programs = identify_stress_programs(loadings, tables.stress_gene_sets)
    stress = stress_scores(perturb_programs, stress_programs)

    perturb_rows = perturb_matrix.reset_index()
    program_rows = perturb_programs.reset_index()
    stress_frame = perturb_rows[["perturbagen", "cell", "time_h", "dose", "tas"]].copy()
    stress_frame["stress_score"] = stress.to_numpy()
    stress_frame["quality_tier"] = stress_frame["tas"].map(lambda x: quality_tier(float(x), config.quality_thresholds))
    stress_frame["stress_z"] = zscore_within_strata(stress_frame, ["cell", "time_h", "quality_tier"], "stress_score")
    stress_frame["stress_high"] = stress_frame["stress_z"] >= config.stress_high_z

    records = []
    perturb_gene_values = perturb_matrix.reset_index()
    for (disease, cohort, method), disease_block in _disease_signature_groups(tables.disease_signatures):
        disease_vec = signature_vector(disease_block, genes)
        up_genes, down_genes = topn_sets(disease_block[disease_block["gene"].isin(genes)], config.top_n)
        disease_program = basis.transform_vector(disease_vec)
        for idx, perturb_row in perturb_gene_values.iterrows():
            metadata = {key: perturb_row[key] for key in ["perturbagen", "cell", "time_h", "dose", "tas"]}
            perturb_vec = perturb_row[genes].to_numpy(dtype=float)
            perturb_series = pd.Series(perturb_vec, index=genes)
            program_vec = program_rows.loc[idx, [c for c in program_rows.columns if c.startswith("program_")]].to_numpy(float)
            tier = quality_tier(float(metadata["tas"]), config.quality_thresholds)
            base = {
                "disease": disease,
                "cohort": cohort,
                "method": method,
                "perturbagen": metadata["perturbagen"],
                "cell": metadata["cell"],
                "time_h": int(metadata["time_h"]),
                "dose": metadata["dose"],
                "tas": float(metadata["tas"]),
                "quality_tier": tier,
            }
            records.append({**base, "objective": "vector", "score": vector_reversal(disease_vec, perturb_vec)})
            records.append({**base, "objective": "connectivity", "score": connectivity_reversal(perturb_series, up_genes, down_genes)})
            records.append({**base, "objective": "program", "score": program_reversal(disease_program, program_vec)})
    rankings = pd.DataFrame.from_records(records)
    rankings = rankings.merge(
        stress_frame[["perturbagen", "cell", "time_h", "dose", "tas", "stress_score", "stress_z", "stress_high"]],
        on=["perturbagen", "cell", "time_h", "dose", "tas"],
        how="left",
    )
    return rankings, loadings, stress_programs


def cross_cohort_stability(rankings: pd.DataFrame, config: OrbitConfig) -> pd.DataFrame:
    rows = []
    group_cols = ["disease", "method", "objective", "cell", "time_h", "quality_tier"]
    for keys, block in rankings.groupby(group_cols, sort=True):
        cohorts = sorted(block["cohort"].unique())
        for a, b in itertools.combinations(cohorts, 2):
            left = block[block["cohort"] == a]
            right = block[block["cohort"] == b]
            metrics = rank_stability(left, right, config.top_k)
            rows.append(dict(zip(group_cols, keys), cohort_a=a, cohort_b=b, **metrics))
    return pd.DataFrame.from_records(rows)


def _aggregate_context(rankings: pd.DataFrame, cells: list[str]) -> pd.DataFrame:
    return (
        rankings[rankings["cell"].isin(cells)]
        .groupby(["disease", "cohort", "method", "objective", "time_h", "quality_tier", "perturbagen"], as_index=False)
        ["score"]
        .median()
    )


def context_sensitivity(rankings: pd.DataFrame, config: OrbitConfig) -> pd.DataFrame:
    rows = []
    for disease, matched in config.matched_contexts.items():
        matched_rankings = _aggregate_context(rankings[rankings["disease"] == disease], matched)
        core_rankings = _aggregate_context(rankings[rankings["disease"] == disease], config.core_contexts)
        merge_cols = ["disease", "cohort", "method", "objective", "time_h", "quality_tier"]
        for keys, left in matched_rankings.groupby(merge_cols, sort=True):
            selector = np.logical_and.reduce([core_rankings[col] == key for col, key in zip(merge_cols, keys)])
            right = core_rankings.loc[selector]
            if right.empty:
                continue
            metrics = rank_stability(left, right, config.top_k)
            rows.append(dict(zip(merge_cols, keys), matched_cells=";".join(matched), core_cells=";".join(config.core_contexts), **metrics))
    return pd.DataFrame.from_records(rows)


def time_stability(rankings: pd.DataFrame, config: OrbitConfig) -> pd.DataFrame:
    rows = []
    group_cols = ["disease", "cohort", "method", "objective", "cell", "quality_tier"]
    for keys, block in rankings.groupby(group_cols, sort=True):
        for a, b in itertools.combinations(sorted(block["time_h"].unique()), 2):
            left = block[block["time_h"] == a]
            right = block[block["time_h"] == b]
            metrics = rank_stability(left, right, config.top_k)
            rows.append(dict(zip(group_cols, keys), time_a=int(a), time_b=int(b), **metrics))
    return pd.DataFrame.from_records(rows)


def variance_decomposition(rankings: pd.DataFrame) -> pd.DataFrame:
    frame = rankings.copy()
    frame["rank_norm_score"] = frame.groupby(["disease", "cohort", "method", "objective", "cell", "time_h", "quality_tier"])["score"].rank(pct=True)
    y = frame["rank_norm_score"].to_numpy(float)
    total = float(np.sum((y - y.mean()) ** 2))
    factors = ["disease", "cohort", "method", "objective", "cell", "time_h", "quality_tier"]
    rows = []
    for factor in factors:
        means = frame.groupby(factor)["rank_norm_score"].transform("mean").to_numpy(float)
        share = float(np.sum((means - y.mean()) ** 2) / total) if total else np.nan
        rows.append({"factor": factor, "variance_share": share})
    used = sum(r["variance_share"] for r in rows if np.isfinite(r["variance_share"]))
    rows.append({"factor": "residual_or_interaction", "variance_share": max(0.0, 1.0 - used)})
    return pd.DataFrame.from_records(rows)


def stress_coupling(rankings: pd.DataFrame, config: OrbitConfig) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    group_cols = ["disease", "cohort", "method", "objective", "cell", "time_h", "quality_tier"]
    for keys, block in rankings.groupby(group_cols, sort=True):
        if len(block) < 3:
            continue
        rho = stats.spearmanr(block["score"], block["stress_z"], nan_policy="omit").statistic
        top = block.sort_values("score", ascending=False).head(max(config.top_k))
        rows.append(
            dict(
                zip(group_cols, keys),
                spearman_reversal_stress=float(rho),
                stress_high_top_fraction=float(top["stress_high"].mean()),
                n=len(block),
            )
        )
    coupling = pd.DataFrame.from_records(rows)

    adjusted = apply_stress_penalty(rankings, config.stress_lambda)
    mitigation_rows = []
    for keys, block in adjusted.groupby(group_cols, sort=True):
        k = max(config.top_k)
        before = block.sort_values("score", ascending=False).head(k)
        after = block.sort_values("adjusted_score", ascending=False).head(k)
        mitigation_rows.append(
            dict(
                zip(group_cols, keys),
                top_k=k,
                stress_high_before=float(before["stress_high"].mean()),
                stress_high_after=float(after["stress_high"].mean()),
                median_score_before=float(before["score"].median()),
                median_score_after=float(after["score"].median()),
            )
        )
    return coupling, pd.DataFrame.from_records(mitigation_rows)


def run_pipeline(tables: OrbitTables, config: OrbitConfig, out_dir: str | Path) -> dict[str, str]:
    out = Path(out_dir)
    tables_dir = out / "tables"
    figures_dir = out / "figures"
    manifest_dir = out / "manifest"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    rankings, loadings, stress_programs = make_rankings(tables, config)
    objective = cross_cohort_stability(rankings, config)
    context = context_sensitivity(rankings, config)
    time = time_stability(rankings, config)
    variance = variance_decomposition(rankings)
    coupling, mitigation = stress_coupling(rankings, config)

    paths = {
        "rankings": tables_dir / "rankings.csv",
        "program_loadings": tables_dir / "program_loadings.csv",
        "stress_programs": tables_dir / "stress_programs.csv",
        "objective_stability": tables_dir / "objective_stability.csv",
        "context_sensitivity": tables_dir / "context_sensitivity.csv",
        "time_stability": tables_dir / "time_stability.csv",
        "variance_decomposition": tables_dir / "variance_decomposition.csv",
        "stress_coupling": tables_dir / "stress_coupling.csv",
        "stress_mitigation": tables_dir / "stress_mitigation.csv",
    }
    rankings.to_csv(paths["rankings"], index=False)
    loadings.to_csv(paths["program_loadings"], index=False)
    stress_programs.to_csv(paths["stress_programs"], index=False)
    objective.to_csv(paths["objective_stability"], index=False)
    context.to_csv(paths["context_sensitivity"], index=False)
    time.to_csv(paths["time_stability"], index=False)
    variance.to_csv(paths["variance_decomposition"], index=False)
    coupling.to_csv(paths["stress_coupling"], index=False)
    mitigation.to_csv(paths["stress_mitigation"], index=False)

    from .reporting import write_figures

    figure_paths = write_figures(objective, context, variance, coupling, figures_dir)
    manifest = write_manifest(out, manifest_dir)
    ro_crate = write_ro_crate(out, out / "ro-crate-metadata.json")
    summary = {
        "rankings_rows": int(len(rankings)),
        "shared_genes": int(
            len(set(tables.disease_signatures["gene"]).intersection(set(tables.perturbation_signatures["gene"])))
        ),
        "objective_stability_rows": int(len(objective)),
        "context_sensitivity_rows": int(len(context)),
        "time_stability_rows": int(len(time)),
        "manifest": str(manifest),
        "ro_crate": str(ro_crate),
        "figures": {key: str(value) for key, value in figure_paths.items()},
    }
    (out / "run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return {key: str(value) for key, value in paths.items()} | {"manifest": str(manifest), "ro_crate": str(ro_crate)}
