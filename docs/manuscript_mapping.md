# Manuscript Mapping

This file maps the ORBIT manuscript sections to repository components.

| Manuscript section | Repository implementation |
| --- | --- |
| 2.3 Benchmark design and disease selection | `configs/*.yaml`, `disease_metadata.csv`, `simulate.py` demo disease panel |
| 2.4 Context panels and exposure time stratification | `OrbitConfig.matched_contexts`, `OrbitConfig.core_contexts`, `context_sensitivity`, `time_stability` |
| 2.5 Data sources and preprocessing | `data.py`, `data/schemas/*.json`, `docs/data_dictionary.md` |
| 2.6 Disease signature construction | `signatures.py` |
| 2.7 Reversal objectives | `objectives.py` |
| 2.8 Robustness evaluation and statistical analysis | `metrics.py`, `cross_cohort_stability`, `variance_decomposition` |
| 2.9 Stress/toxicity confounding and mitigation | `stress.py`, `stress_coupling` |
| 2.10 Null controls | `metrics.py`, `simulate.py`, and user-supplied randomized/sign-flipped signature tables |
| Results figures and tables | `reporting.py`, `results/*/tables`, `results/*/figures` |
| FAIR availability | `README.md`, `CITATION.cff`, `codemeta.json`, `docs/fairness_checklist.md`, `fair.py` |

The code uses a consistent score orientation: higher scores indicate stronger reversal for all objectives. This makes ranking and top-K reporting uniform across vector, connectivity-style, and program-level objectives.

