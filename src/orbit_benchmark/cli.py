from __future__ import annotations

import argparse
import json
from pathlib import Path

from .config import load_config
from .data import load_tables, validate_input_dir
from .evaluation import run_pipeline
from .fair import write_manifest, write_ro_crate
from .simulate import write_example_dataset


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="orbit", description="Run the ORBIT transcriptome inversion benchmark.")
    sub = parser.add_subparsers(dest="command", required=True)

    init = sub.add_parser("init-example", help="Create deterministic example input tables.")
    init.add_argument("--out", required=True, help="Output directory for CSV inputs.")
    init.add_argument("--seed", type=int, default=2026)
    init.add_argument("--genes", type=int, default=260)
    init.add_argument("--perturbagens", type=int, default=96)

    validate = sub.add_parser("validate", help="Validate an ORBIT input directory.")
    validate.add_argument("--data", required=True)

    run = sub.add_parser("run", help="Run the benchmark pipeline.")
    run.add_argument("--data", required=True)
    run.add_argument("--config", default=None)
    run.add_argument("--out", required=True)

    manifest = sub.add_parser("manifest", help="Create file checksum manifest and RO-Crate metadata.")
    manifest.add_argument("--root", required=True)
    manifest.add_argument("--out", required=True)

    fair = sub.add_parser("fair-check", help="Validate inputs and write a manifest for inputs and outputs.")
    fair.add_argument("--data", required=True)
    fair.add_argument("--out", required=True)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "init-example":
        write_example_dataset(args.out, seed=args.seed, n_genes=args.genes, n_perturbagens=args.perturbagens)
        print(f"example data written to {Path(args.out).resolve()}")
        return 0

    if args.command == "validate":
        summary = validate_input_dir(args.data)
        print(json.dumps(summary, indent=2))
        return 0

    if args.command == "run":
        config = load_config(args.config)
        tables = load_tables(args.data)
        outputs = run_pipeline(tables, config, args.out)
        print(json.dumps(outputs, indent=2))
        return 0

    if args.command == "manifest":
        manifest_path = write_manifest(args.root, args.out)
        ro_crate = write_ro_crate(args.root, Path(args.out) / "ro-crate-metadata.json")
        print(json.dumps({"manifest": str(manifest_path), "ro_crate": str(ro_crate)}, indent=2))
        return 0

    if args.command == "fair-check":
        summary = validate_input_dir(args.data)
        manifest_path = write_manifest(args.out, Path(args.out) / "manifest")
        data_manifest = write_manifest(args.data, Path(args.out) / "input_manifest")
        payload = {"validation": summary, "output_manifest": str(manifest_path), "input_manifest": str(data_manifest)}
        print(json.dumps(payload, indent=2))
        return 0

    parser.error(f"unknown command: {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())

