from pathlib import Path

import pandas as pd

from orbit_benchmark.config import OrbitConfig
from orbit_benchmark.data import load_tables, validate_input_dir
from orbit_benchmark.evaluation import run_pipeline
from orbit_benchmark.simulate import write_example_dataset


def test_demo_pipeline_runs(tmp_path: Path):
    data_dir = tmp_path / "data"
    out_dir = tmp_path / "results"
    write_example_dataset(data_dir, seed=7, n_genes=80, n_perturbagens=24)
    summary = validate_input_dir(data_dir)
    assert summary["diseases"] == 3
    config = OrbitConfig(
        top_n=12,
        top_k=(5, 10),
        program_k=6,
        matched_contexts={
            "psoriasis": ["A375", "HA1E"],
            "er_breast_cancer": ["MCF7", "BT20"],
            "lung_adenocarcinoma": ["A549", "PC3"],
        },
        core_contexts=["A375", "A549", "BT20", "HA1E", "HT29", "MCF7", "PC3"],
        bootstrap_iterations=5,
    )
    run_pipeline(load_tables(data_dir), config, out_dir)
    objective = pd.read_csv(out_dir / "tables" / "objective_stability.csv")
    assert {"objective", "kendall_tau", "jaccard_at_5"}.issubset(objective.columns)
    assert (out_dir / "manifest" / "files_manifest.csv").exists()

