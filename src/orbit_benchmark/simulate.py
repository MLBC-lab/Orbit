from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .signatures import benjamini_hochberg


DISEASES = ["psoriasis", "er_breast_cancer", "lung_adenocarcinoma"]
COHORTS = ["cohort_a", "cohort_b", "cohort_c"]
CELLS = ["A375", "A549", "BT20", "HA1E", "HT29", "MCF7", "PC3"]
TIMES = [6, 24]


def _two_sided_p_from_z(z: np.ndarray) -> np.ndarray:
    from scipy import stats

    return 2.0 * stats.norm.sf(np.abs(z))


def write_example_dataset(out_dir: str | Path, *, seed: int = 2026, n_genes: int = 260, n_perturbagens: int = 96) -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    genes = np.array([f"GENE{i:04d}" for i in range(1, n_genes + 1)])
    n_programs = 16
    loadings = rng.normal(0, 0.09, size=(n_programs, n_genes))
    stress_gene_idx = np.unique(
        np.r_[
            0 : min(20, n_genes),
            min(40, n_genes) : min(50, n_genes),
            min(75, n_genes) : min(85, n_genes),
        ]
    )
    if len(stress_gene_idx) < min(8, n_genes):
        stress_gene_idx = np.arange(min(8, n_genes))
    loadings[0, stress_gene_idx] += rng.normal(0.55, 0.05, size=len(stress_gene_idx))
    loadings[1, stress_gene_idx[:20]] -= rng.normal(0.42, 0.04, size=20)

    disease_programs = {
        "psoriasis": np.array([0.2, 0.1, 1.2, -0.5, 0.7, 0.0, -0.2, 0.3, 0, 0, 0, 0, 0, 0, 0, 0]),
        "er_breast_cancer": np.array([0.1, -0.1, -0.4, 1.4, 0.2, 0.6, 0.1, -0.3, 0, 0, 0, 0, 0, 0, 0, 0]),
        "lung_adenocarcinoma": np.array([0.15, 0.2, 0.2, -0.2, 1.1, -0.5, 0.8, 0.4, 0, 0, 0, 0, 0, 0, 0, 0]),
    }

    disease_rows = []
    for disease in DISEASES:
        base = disease_programs[disease] @ loadings
        for cohort in COHORTS:
            noise = rng.normal(0, 0.28, size=n_genes)
            effect = base + noise
            p_value = _two_sided_p_from_z(effect / 0.35)
            q_value = benjamini_hochberg(p_value)
            for gene, eff, p, q in zip(genes, effect, p_value, q_value):
                disease_rows.append(
                    {
                        "disease": disease,
                        "cohort": cohort,
                        "method": "cohort_de",
                        "gene": gene,
                        "effect": float(eff),
                        "p_value": float(p),
                        "q_value": float(q),
                    }
                )

    perturbagens = [f"ORB{idx:04d}" for idx in range(1, n_perturbagens + 1)]
    perturb_meta = []
    perturb_rows = []
    for i, perturbagen in enumerate(perturbagens):
        family = ["candidate", "stress_like", "weak", "disease_mimic"][i % 4]
        perturb_meta.append({"perturbagen": perturbagen, "family": family, "touchstone": True})
        for cell in CELLS:
            cell_shift = rng.normal(0, 0.18, size=n_programs)
            cell_shift[0] += {"A375": 0.2, "A549": 0.1, "MCF7": -0.1}.get(cell, 0.0)
            for time_h in TIMES:
                time_shift = rng.normal(0, 0.12, size=n_programs)
                time_shift[0] += 0.25 if time_h == 24 else -0.05
                activ = rng.normal(0, 0.8, size=n_programs) + cell_shift + time_shift
                if family == "candidate":
                    disease = DISEASES[(i // 4) % len(DISEASES)]
                    activ -= 0.65 * disease_programs[disease]
                elif family == "stress_like":
                    activ[0] += 2.1
                    activ[1] -= 0.8
                elif family == "disease_mimic":
                    disease = DISEASES[(i // 4) % len(DISEASES)]
                    activ += 0.45 * disease_programs[disease]
                tas = float(np.clip(0.08 + 0.18 * np.linalg.norm(activ[:8]) + rng.normal(0, 0.09), 0.02, 0.95))
                z = activ @ loadings + rng.normal(0, 0.22, size=n_genes)
                for gene, value in zip(genes, z):
                    perturb_rows.append(
                        {
                            "perturbagen": perturbagen,
                            "cell": cell,
                            "time_h": time_h,
                            "dose": "10uM",
                            "tas": tas,
                            "gene": gene,
                            "z": float(value),
                        }
                    )

    stress_sets = []
    stress_genes = genes[stress_gene_idx]
    for set_name, members in {
        "hallmark_apoptosis_like": stress_genes[:24],
        "hallmark_tnfa_nfkb_like": stress_genes[10:34],
        "hallmark_unfolded_protein_response_like": stress_genes[20:],
    }.items():
        for gene in members:
            stress_sets.append({"set_name": set_name, "gene": str(gene)})

    disease_meta = pd.DataFrame(
        {
            "disease": DISEASES,
            "category": ["inflammatory", "solid_tumor", "solid_tumor"],
            "selection_note": [
                "multi-cohort inflammatory benchmark",
                "hormone-positive cancer benchmark",
                "epithelial cancer benchmark",
            ],
        }
    )

    pd.DataFrame(disease_rows).to_csv(out / "disease_signatures.csv", index=False)
    pd.DataFrame(perturb_rows).to_csv(out / "perturbation_signatures.csv", index=False)
    pd.DataFrame(stress_sets).to_csv(out / "stress_gene_sets.csv", index=False)
    disease_meta.to_csv(out / "disease_metadata.csv", index=False)
    pd.DataFrame(perturb_meta).to_csv(out / "perturbation_metadata.csv", index=False)
