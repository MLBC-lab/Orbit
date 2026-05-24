from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def _save(fig: plt.Figure, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def write_figures(
    objective_stability: pd.DataFrame,
    context_sensitivity: pd.DataFrame,
    variance_decomposition: pd.DataFrame,
    stress_coupling: pd.DataFrame,
    out_dir: str | Path,
) -> dict[str, Path]:
    out = Path(out_dir)
    paths: dict[str, Path] = {}

    fig, ax = plt.subplots(figsize=(7, 4.5))
    if not objective_stability.empty:
        jaccard_cols = sorted(
            [c for c in objective_stability.columns if c.startswith("jaccard_at_")],
            key=lambda c: int(c.rsplit("_", 1)[1]),
        )
        jaccard_col = jaccard_cols[-1] if jaccard_cols else "kendall_tau"
        summary = objective_stability.groupby("objective", as_index=False).agg(
            kendall_tau=("kendall_tau", "median"),
            topk_jaccard=(jaccard_col, "median"),
        )
        for _, row in summary.iterrows():
            ax.scatter(row["kendall_tau"], row["topk_jaccard"], s=90, label=row["objective"])
        ax.legend(frameon=False)
    ax.set_xlabel("Median Kendall tau")
    ax.set_ylabel("Median top-K Jaccard")
    ax.set_title("Objective stability")
    ax.grid(alpha=0.25)
    paths["objective_stability"] = _save(fig, out / "objective_stability.png")

    fig, ax = plt.subplots(figsize=(7, 4.5))
    if not context_sensitivity.empty:
        summary = context_sensitivity.groupby("objective", as_index=False)["kendall_tau"].median()
        ax.bar(summary["objective"], summary["kendall_tau"], color=["#5B8FF9", "#61DDAA", "#F6BD16"][: len(summary)])
    ax.set_ylabel("Matched vs core Kendall tau")
    ax.set_title("Context sensitivity")
    ax.grid(axis="y", alpha=0.25)
    paths["context_sensitivity"] = _save(fig, out / "context_sensitivity.png")

    fig, ax = plt.subplots(figsize=(7, 4.5))
    if not variance_decomposition.empty:
        plot_data = variance_decomposition.sort_values("variance_share", ascending=False)
        ax.barh(plot_data["factor"], plot_data["variance_share"], color="#5AD8A6")
        ax.invert_yaxis()
    ax.set_xlabel("Variance share")
    ax.set_title("Benchmark factor variance shares")
    ax.grid(axis="x", alpha=0.25)
    paths["variance_decomposition"] = _save(fig, out / "variance_decomposition.png")

    fig, ax = plt.subplots(figsize=(7, 4.5))
    if not stress_coupling.empty:
        summary = stress_coupling.groupby("objective", as_index=False)["spearman_reversal_stress"].median()
        ax.bar(summary["objective"], summary["spearman_reversal_stress"], color=["#9270CA", "#FF9D4D", "#269A99"][: len(summary)])
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("Median Spearman rho")
    ax.set_title("Reversal and stress coupling")
    ax.grid(axis="y", alpha=0.25)
    paths["stress_coupling"] = _save(fig, out / "stress_coupling.png")

    return paths
