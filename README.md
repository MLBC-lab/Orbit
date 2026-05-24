# ORBIT Benchmark

ORBIT is a reference implementation of the **Objective Robustness Benchmark for Inversion of Transcriptomes**. It evaluates how disease-signature reversal rankings change across disease cohorts, reversal objectives, cellular context, exposure time, perturbation quality, and generic stress activation.

The implementation follows the manuscript structure:

- disease signature construction and consensus summaries
- vector, connectivity-style, and program-level reversal objectives
- cross-cohort, cross-context, and 6 h vs 24 h rank stability
- top-K agreement with Jaccard and overlap coefficient
- variance-share summaries for benchmark factors
- stress-confounding diagnostics and stress-penalized reranking
- prospective null controls and bootstrap uncertainty
- FAIR metadata, manifests, schemas, and checksum generation

The repository does not bundle restricted CMap/LINCS matrices. A deterministic ORBIT-like example dataset is generated locally so that installation, scoring, reporting, and tests can be run immediately.

## Quick Start

```powershell
cd C:\Writing\School\Orbit\orbit-benchmark
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
orbit init-example --out data/example --seed 2026
orbit run --data data/example --config configs/demo.yaml --out results/demo
orbit fair-check --data data/example --out results/demo
pytest
```

On macOS or Linux, replace the activation line with:

```bash
source .venv/bin/activate
```

The main outputs are written to `results/demo`:

- `tables/objective_stability.csv`
- `tables/context_sensitivity.csv`
- `tables/time_stability.csv`
- `tables/variance_decomposition.csv`
- `tables/stress_coupling.csv`
- `tables/stress_mitigation.csv`
- `figures/objective_stability.png`
- `figures/context_sensitivity.png`
- `figures/variance_decomposition.png`
- `figures/stress_coupling.png`
- `manifest/files_manifest.csv`
- `run_summary.json`

## Input Tables

Real analyses use the same table layout as the demo data. All files are plain CSV with UTF-8 encoding.

### `disease_signatures.csv`

One row per disease, cohort, method, and gene.

| column | description |
| --- | --- |
| `disease` | Disease label. |
| `cohort` | Independent disease cohort or cohort-consensus label. |
| `method` | Signature method, such as `limma_topn` or `stouffer_consensus`. |
| `gene` | HGNC symbol or stable feature identifier. |
| `effect` | Signed disease effect statistic. Positive values are up in disease. |
| `p_value` | Nominal differential-expression p-value. |
| `q_value` | Adjusted p-value. |

### `perturbation_signatures.csv`

One row per perturbation signature and gene.

| column | description |
| --- | --- |
| `perturbagen` | Compound, genetic perturbation, or perturbagen identifier. |
| `cell` | Cell context. |
| `time_h` | Exposure time in hours. |
| `dose` | Dose label. |
| `tas` | Transcriptional Activity Score. |
| `gene` | Landmark gene or aligned feature identifier. |
| `z` | Perturbational signature value. |

### `stress_gene_sets.csv`

One row per stress-related gene-set membership.

| column | description |
| --- | --- |
| `set_name` | Stress-related pathway or collection label. |
| `gene` | Gene identifier in the same namespace as the signatures. |

### Optional Metadata

`disease_metadata.csv` and `perturbation_metadata.csv` can be supplied to make reports more interpretable. The pipeline will run without them.

## Configuration

`configs/demo.yaml` is tuned for fast local execution. The more paper-like settings are in `configs/paper_like.yaml`.

Important fields:

- `top_n`: genes selected per direction for enrichment objectives
- `top_k`: top-K stability cutoffs
- `program_k`: number of PCA programs
- `quality_thresholds`: TAS strata
- `matched_contexts`: disease-specific matched cell panels
- `core_contexts`: larger cell panel used for context sensitivity
- `stress_lambda`: penalty used for stress-aware reranking
- `bootstrap_iterations`: number of bootstrap replicates for confidence intervals

## Running with Real ORBIT Data

1. Export or construct disease signatures into `disease_signatures.csv`.
2. Export CMap/LINCS Level-5 landmark signatures into `perturbation_signatures.csv`.
3. Keep perturbagens in comparable dose, cell, time, and TAS strata before ranking.
4. Add stress-related gene sets to `stress_gene_sets.csv`.
5. Run:

```powershell
orbit validate --data path\to\orbit_tables
orbit run --data path\to\orbit_tables --config configs/paper_like.yaml --out results\paper_like
orbit manifest --root results\paper_like --out results\paper_like\manifest
```

## FAIR Practices

This repository includes:

- `CITATION.cff` for citation metadata
- `codemeta.json` for machine-readable software metadata
- `LICENSE` for reuse terms
- `data/schemas/*.json` for table interoperability
- `docs/data_dictionary.md` for semantic definitions
- `docs/fairness_checklist.md` for FAIR review
- `orbit manifest` for checksums and file-level provenance
- deterministic demo generation through explicit seeds

The code is intended to be run as written. Any analysis using non-public or restricted matrices should provide accession identifiers, preprocessing notes, and checksum manifests for all derived tables.

## Development Checks

```powershell
pytest
orbit init-example --out data/example --seed 2026
orbit validate --data data/example
orbit run --data data/example --config configs/demo.yaml --out results/demo
```

## Repository Layout

```text
configs/              Pipeline configurations
data/schemas/         JSON schemas for required CSV inputs
docs/                 Data dictionary and FAIR checklist
scripts/              Thin wrappers for common workflows
src/orbit_benchmark/  Python package
tests/                Unit and pipeline smoke tests
```

