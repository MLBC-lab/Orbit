from __future__ import annotations

from pathlib import Path

from orbit_benchmark.config import load_config
from orbit_benchmark.data import load_tables
from orbit_benchmark.evaluation import run_pipeline
from orbit_benchmark.simulate import write_example_dataset


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    data_dir = root / "data" / "example"
    out_dir = root / "results" / "demo"
    write_example_dataset(data_dir, seed=2026)
    config = load_config(root / "configs" / "demo.yaml")
    run_pipeline(load_tables(data_dir), config, out_dir)
    print(f"demo results written to {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

