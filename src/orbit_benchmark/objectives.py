from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def cosine_similarity(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    if denom == 0:
        return 0.0
    return float(np.dot(a, b) / denom)


def vector_reversal(disease: np.ndarray, perturbation: np.ndarray) -> float:
    return -cosine_similarity(disease, perturbation)


def weighted_enrichment_score(ranked_genes: list[str], gene_set: set[str], weights: np.ndarray | None = None) -> float:
    if not gene_set:
        return 0.0
    genes = np.asarray(ranked_genes, dtype=object)
    if weights is None:
        weights = np.ones(len(genes), dtype=float)
    else:
        weights = np.abs(np.asarray(weights, dtype=float))
    hits = np.isin(genes, list(gene_set))
    hit_count = int(hits.sum())
    if hit_count == 0 or hit_count == len(genes):
        return 0.0
    hit_weights = weights * hits
    hit_norm = hit_weights.sum()
    miss_norm = len(genes) - hit_count
    running = np.cumsum(np.where(hits, hit_weights / hit_norm, -1.0 / miss_norm))
    idx = int(np.argmax(np.abs(running)))
    return float(running[idx])


def connectivity_reversal(perturbation: pd.Series, up_genes: set[str], down_genes: set[str]) -> float:
    ranked = perturbation.sort_values(ascending=False)
    genes = [str(g) for g in ranked.index]
    weights = ranked.to_numpy(dtype=float)
    es_up = weighted_enrichment_score(genes, up_genes, weights)
    es_down = weighted_enrichment_score(genes, down_genes, weights)
    return float((es_down - es_up) / 2.0)


@dataclass(frozen=True)
class ProgramBasis:
    genes: list[str]
    scaler: StandardScaler
    pca: PCA

    def transform_vector(self, vector: np.ndarray) -> np.ndarray:
        frame = pd.DataFrame([vector], columns=self.genes)
        return self.pca.transform(self.scaler.transform(frame))[0]

    def transform_matrix(self, matrix: pd.DataFrame) -> pd.DataFrame:
        aligned = matrix.reindex(columns=self.genes).fillna(0.0)
        values = self.pca.transform(self.scaler.transform(aligned))
        columns = [f"program_{i + 1:03d}" for i in range(values.shape[1])]
        return pd.DataFrame(values, index=aligned.index, columns=columns)

    def loadings(self) -> pd.DataFrame:
        rows = []
        for i, component in enumerate(self.pca.components_, start=1):
            for gene, loading in zip(self.genes, component):
                rows.append({"program": f"program_{i:03d}", "gene": gene, "loading": float(loading)})
        return pd.DataFrame(rows)


def learn_program_basis(matrix: pd.DataFrame, k: int, random_seed: int = 2026) -> ProgramBasis:
    genes = [str(g) for g in matrix.columns]
    n_components = max(1, min(k, matrix.shape[0] - 1, matrix.shape[1]))
    scaler = StandardScaler()
    scaled = scaler.fit_transform(matrix.reindex(columns=genes).fillna(0.0))
    pca = PCA(n_components=n_components, random_state=random_seed)
    pca.fit(scaled)
    return ProgramBasis(genes=genes, scaler=scaler, pca=pca)


def program_reversal(disease_programs: np.ndarray, perturbation_programs: np.ndarray) -> float:
    return -cosine_similarity(disease_programs, perturbation_programs)

