from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd


DISEASE_COLUMNS = {"disease", "cohort", "method", "gene", "effect", "p_value", "q_value"}
PERTURBATION_COLUMNS = {"perturbagen", "cell", "time_h", "dose", "tas", "gene", "z"}
STRESS_COLUMNS = {"set_name", "gene"}


@dataclass(frozen=True)
class OrbitTables:
    disease_signatures: pd.DataFrame
    perturbation_signatures: pd.DataFrame
    stress_gene_sets: pd.DataFrame
    disease_metadata: pd.DataFrame | None = None
    perturbation_metadata: pd.DataFrame | None = None


def read_table(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path)


def validate_columns(frame: pd.DataFrame, required: set[str], name: str) -> None:
    missing = sorted(required.difference(frame.columns))
    if missing:
        joined = ", ".join(missing)
        raise ValueError(f"{name} is missing required columns: {joined}")


def validate_input_dir(data_dir: str | Path) -> dict[str, int]:
    root = Path(data_dir)
    disease = read_table(root / "disease_signatures.csv")
    perturb = read_table(root / "perturbation_signatures.csv")
    stress = read_table(root / "stress_gene_sets.csv")
    validate_columns(disease, DISEASE_COLUMNS, "disease_signatures.csv")
    validate_columns(perturb, PERTURBATION_COLUMNS, "perturbation_signatures.csv")
    validate_columns(stress, STRESS_COLUMNS, "stress_gene_sets.csv")
    return {
        "disease_signature_rows": int(len(disease)),
        "perturbation_signature_rows": int(len(perturb)),
        "stress_gene_set_rows": int(len(stress)),
        "diseases": int(disease["disease"].nunique()),
        "perturbagens": int(perturb["perturbagen"].nunique()),
        "genes": int(len(set(disease["gene"]).intersection(set(perturb["gene"])))),
    }


def load_tables(data_dir: str | Path) -> OrbitTables:
    root = Path(data_dir)
    disease = read_table(root / "disease_signatures.csv")
    perturb = read_table(root / "perturbation_signatures.csv")
    stress = read_table(root / "stress_gene_sets.csv")
    validate_columns(disease, DISEASE_COLUMNS, "disease_signatures.csv")
    validate_columns(perturb, PERTURBATION_COLUMNS, "perturbation_signatures.csv")
    validate_columns(stress, STRESS_COLUMNS, "stress_gene_sets.csv")

    disease_meta = None
    if (root / "disease_metadata.csv").exists():
        disease_meta = read_table(root / "disease_metadata.csv")

    perturb_meta = None
    if (root / "perturbation_metadata.csv").exists():
        perturb_meta = read_table(root / "perturbation_metadata.csv")

    disease["gene"] = disease["gene"].astype(str)
    perturb["gene"] = perturb["gene"].astype(str)
    stress["gene"] = stress["gene"].astype(str)
    perturb["time_h"] = perturb["time_h"].astype(int)
    return OrbitTables(disease, perturb, stress, disease_meta, perturb_meta)


def pivot_perturbations(perturbation_signatures: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    meta_cols = ["perturbagen", "cell", "time_h", "dose", "tas"]
    meta = perturbation_signatures[meta_cols].drop_duplicates().reset_index(drop=True)
    table = perturbation_signatures.pivot_table(
        index=meta_cols,
        columns="gene",
        values="z",
        aggfunc="mean",
    )
    table = table.sort_index(axis=1).fillna(0.0)
    return table, meta

