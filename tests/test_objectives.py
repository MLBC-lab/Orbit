import numpy as np
import pandas as pd

from orbit_benchmark.objectives import connectivity_reversal, vector_reversal, weighted_enrichment_score


def test_vector_reversal_prefers_opposite_direction():
    disease = np.array([1.0, 0.5, -1.0])
    opposite = np.array([-1.0, -0.5, 1.0])
    same = np.array([1.0, 0.5, -1.0])
    assert vector_reversal(disease, opposite) > vector_reversal(disease, same)


def test_weighted_enrichment_score_sign_tracks_rank_position():
    genes = ["A", "B", "C", "D", "E"]
    weights = np.ones(5)
    assert weighted_enrichment_score(genes, {"A", "B"}, weights) > 0
    assert weighted_enrichment_score(genes, {"D", "E"}, weights) < 0


def test_connectivity_reversal_rewards_down_up_inversion():
    perturb = pd.Series({"UP1": -2.0, "UP2": -1.8, "MID": 0.0, "DN1": 1.7, "DN2": 2.0})
    score = connectivity_reversal(perturb, {"UP1", "UP2"}, {"DN1", "DN2"})
    assert score > 0

