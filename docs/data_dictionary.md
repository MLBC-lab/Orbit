# Data Dictionary

## Disease Signature

`effect` is a signed statistic for case-control differential expression. Positive values indicate genes up-regulated in disease relative to control. Negative values indicate genes down-regulated in disease.

`method` records how the disease signature was constructed. Common values include:

- `cohort_de`: cohort-level differential expression
- `limma_topn`: top-N signature from a moderated linear model
- `characteristic_direction`: characteristic-direction sensitivity signature
- `stouffer_consensus`: consensus signature from signed Z-score aggregation
- `rra_consensus`: robust rank aggregation consensus

## Perturbation Signature

`z` is the signed perturbational response for one gene in a specific perturbagen, cell context, exposure time, dose, and quality stratum.

`tas` is the Transcriptional Activity Score used to stratify signature quality:

- low or noisy: `tas < 0.2`
- moderate: `tas >= 0.2`
- stringent: `tas >= 0.5`

## Reversal Score Direction

This implementation reports higher values as stronger disease reversal for every objective. The connectivity-style score is therefore oriented as `(ES_down - ES_up) / 2`, where `ES_down` measures disease-down genes enriched at the top of the perturbational ranking and `ES_up` measures disease-up genes enriched at the top.

## Robustness Metrics

`kendall_tau` measures pairwise ordering agreement between two rankings.

`spearman_rho` measures rank correlation using rank differences.

`jaccard_at_k` measures overlap among the top-K perturbagens.

`overlap_at_k` divides the top-K intersection by the smaller top-K set size.

## Stress Diagnostics

Stress-like programs are identified by enrichment of stress-related gene sets among high-magnitude PCA loadings. Perturbagens are labeled stress-high when the within-stratum stress score is at least `stress_high_z` standard deviations above the stratum mean.

