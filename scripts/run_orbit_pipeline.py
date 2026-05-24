from __future__ import annotations

import argparse

from orbit_benchmark.config import load_config
from orbit_benchmark.data import load_tables
from orbit_benchmark.evaluation import run_pipeline


def main() -> int:
    parser = argparse.ArgumentParser(description="Run ORBIT against prepared CSV input tables.")
    parser.add_argument("--data", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    run_pipeline(load_tables(args.data), load_config(args.config), args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

