import pandas as pd

from orbit_benchmark.metrics import jaccard_at_k, rank_stability


def test_topk_overlap_uses_perturbagen_identity():
    a = pd.DataFrame({"perturbagen": ["a", "b", "c"], "score": [3, 2, 1]})
    b = pd.DataFrame({"perturbagen": ["b", "a", "c"], "score": [3, 2, 1]})
    assert jaccard_at_k(a, b, 2) == 1.0


def test_rank_stability_reports_common_count():
    a = pd.DataFrame({"perturbagen": ["a", "b", "c"], "score": [3, 2, 1]})
    b = pd.DataFrame({"perturbagen": ["a", "b", "c"], "score": [3, 2, 1]})
    metrics = rank_stability(a, b, (2,))
    assert metrics["n_common"] == 3
    assert metrics["kendall_tau"] > 0.99

