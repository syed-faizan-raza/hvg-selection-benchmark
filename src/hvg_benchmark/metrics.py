from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import LabelEncoder

from .hvg import select_top_genes


def _encoded(labels) -> np.ndarray:
    return LabelEncoder().fit_transform(pd.Series(labels).astype(str))


def prepare_embedding(
    adata,
    genes: Iterable[str],
    *,
    counts_layer: str = "counts",
    n_pcs: int = 30,
    n_neighbors: int = 15,
    random_state: int = 7,
    compute_umap: bool = False,
):
    """Normalize, log-transform, PCA, and build the neighbor graph on selected genes."""
    genes = [gene for gene in genes if gene in adata.var_names]
    if len(genes) < 2:
        raise ValueError("At least two selected genes are required to compute an embedding.")

    ad = adata[:, genes].copy()
    if counts_layer in ad.layers:
        ad.X = ad.layers[counts_layer].copy()

    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)

    n_comps = min(int(n_pcs), ad.n_obs - 1, ad.n_vars - 1)
    n_comps = max(2, n_comps)
    sc.tl.pca(ad, n_comps=n_comps, svd_solver="arpack", random_state=random_state)
    sc.pp.neighbors(
        ad,
        n_neighbors=min(int(n_neighbors), ad.n_obs - 1),
        n_pcs=n_comps,
        random_state=random_state,
    )
    if compute_umap:
        sc.tl.umap(ad, random_state=random_state)
    return ad


def sweep_leiden_best(
    adata,
    labels,
    *,
    resolutions: Iterable[float] = tuple(np.round(np.arange(0.1, 2.1, 0.1), 2)),
    random_state: int = 7,
) -> dict[str, float]:
    """Sweep Leiden resolution and keep the best ARI and NMI against labels."""
    y_true = _encoded(labels)
    best = {
        "ARI": np.nan,
        "NMI": np.nan,
        "best_ari_resolution": np.nan,
        "best_nmi_resolution": np.nan,
    }

    for resolution in resolutions:
        key = "leiden_tmp"
        sc.tl.leiden(
            adata,
            resolution=float(resolution),
            key_added=key,
            random_state=random_state,
        )
        y_pred = _encoded(adata.obs[key])
        ari = adjusted_rand_score(y_true, y_pred)
        nmi = normalized_mutual_info_score(y_true, y_pred)
        if np.isnan(best["ARI"]) or ari > best["ARI"]:
            best["ARI"] = float(ari)
            best["best_ari_resolution"] = float(resolution)
        if np.isnan(best["NMI"]) or nmi > best["NMI"]:
            best["NMI"] = float(nmi)
            best["best_nmi_resolution"] = float(resolution)
    return best


def knn_accuracy(adata, labels, *, k: int = 5) -> float:
    """Majority-label kNN accuracy in PCA space."""
    X = adata.obsm["X_pca"]
    y = _encoded(labels)
    k = min(int(k), X.shape[0] - 1)
    nbrs = NearestNeighbors(n_neighbors=k + 1).fit(X)
    indices = nbrs.kneighbors(X, return_distance=False)[:, 1:]

    correct = 0
    for i, neigh in enumerate(indices):
        values, counts = np.unique(y[neigh], return_counts=True)
        correct += int(values[np.argmax(counts)] == y[i])
    return float(correct / len(y))


def average_silhouette_width(adata, labels, *, max_cells: int = 5000, random_state: int = 7) -> float:
    """ASW on PCA coordinates, sampled for large datasets."""
    X = adata.obsm["X_pca"]
    y = _encoded(labels)
    if len(np.unique(y)) < 2:
        return np.nan
    sample_size = min(max_cells, X.shape[0]) if X.shape[0] > max_cells else None
    return float(
        silhouette_score(X, y, metric="euclidean", sample_size=sample_size, random_state=random_state)
    )


def variance_ratio(adata, labels) -> float:
    """Between-label variance divided by within-label variance in PCA space."""
    X = np.asarray(adata.obsm["X_pca"])
    y = _encoded(labels)
    grand = X.mean(axis=0)

    between = 0.0
    within = 0.0
    for label in np.unique(y):
        Xg = X[y == label]
        center = Xg.mean(axis=0)
        between += Xg.shape[0] * float(np.sum((center - grand) ** 2))
        within += float(np.sum((Xg - center) ** 2))
    return float(between / (within + 1e-12))


def _lisi_components(adata, labels, *, k: int = 30) -> tuple[np.ndarray, np.ndarray]:
    X = np.asarray(adata.obsm["X_pca"])
    y = _encoded(labels)
    k = min(int(k), X.shape[0] - 1)
    indices = NearestNeighbors(n_neighbors=k + 1).fit(X).kneighbors(X, return_distance=False)[:, 1:]

    inverse_values = []
    lisi_values = []
    for neigh in indices:
        _, counts = np.unique(y[neigh], return_counts=True)
        p = counts / counts.sum()
        inverse = float(np.sum(p**2))
        inverse_values.append(inverse)
        lisi_values.append(float(1.0 / inverse))
    return np.asarray(inverse_values), np.asarray(lisi_values)


def inverse_lisi(adata, labels, *, k: int = 30) -> float:
    """Mean inverse LISI, where larger values indicate purer local neighborhoods."""
    inverse_values, _ = _lisi_components(adata, labels, k=k)
    return float(np.mean(inverse_values))


def lisi_score(adata, labels, *, k: int = 30) -> float:
    """Mean LISI, where larger values indicate stronger local mixing."""
    _, lisi_values = _lisi_components(adata, labels, k=k)
    return float(np.mean(lisi_values))


def evaluate_method_at_k(
    adata,
    scores: pd.Series,
    *,
    dataset: str,
    method: str,
    k: int,
    label_key: str | None,
    batch_key: str | None = None,
    counts_layer: str = "counts",
    n_pcs: int = 30,
    n_neighbors: int = 15,
    random_state: int = 7,
):
    """Evaluate one method at one HVG count."""
    genes = select_top_genes(scores, k)
    ad_emb = prepare_embedding(
        adata,
        genes,
        counts_layer=counts_layer,
        n_pcs=n_pcs,
        n_neighbors=n_neighbors,
        random_state=random_state,
    )

    row = {
        "dataset": dataset,
        "method": method,
        "K": int(k),
        "n_cells": int(ad_emb.n_obs),
        "n_genes_used": int(ad_emb.n_vars),
        "n_pcs": int(ad_emb.obsm["X_pca"].shape[1]),
    }

    if batch_key is not None and batch_key in adata.obs:
        batches = adata.obs[batch_key].values
        row["batch_LISI"] = lisi_score(ad_emb, batches)
    else:
        row["batch_LISI"] = np.nan

    if label_key is None or label_key not in adata.obs:
        row.update(
            ARI=np.nan,
            NMI=np.nan,
            best_ari_resolution=np.nan,
            best_nmi_resolution=np.nan,
            kNN_acc=np.nan,
            ASW=np.nan,
            VarRatio=np.nan,
            invLISI=np.nan,
        )
        return row

    labels = adata.obs[label_key].values
    row.update(sweep_leiden_best(ad_emb, labels, random_state=random_state))
    row.update(
        kNN_acc=knn_accuracy(ad_emb, labels, k=5),
        ASW=average_silhouette_width(ad_emb, labels, random_state=random_state),
        VarRatio=variance_ratio(ad_emb, labels),
        invLISI=inverse_lisi(ad_emb, labels),
    )
    return row


def evaluate_scores_grid(
    adata,
        scores_by_method: dict[str, pd.Series],
    *,
    dataset: str,
    k_values: Iterable[int],
    label_key: str | None,
    batch_key: str | None = None,
    counts_layer: str = "counts",
    n_pcs: int = 30,
    n_neighbors: int = 15,
    random_state: int = 7,
) -> pd.DataFrame:
    """Evaluate every method/K pair and return a tidy metrics table."""
    rows = []
    for method, scores in scores_by_method.items():
        for k in k_values:
            rows.append(
                evaluate_method_at_k(
                    adata,
                    scores,
                    dataset=dataset,
                    method=method,
                    k=int(k),
                    label_key=label_key,
                    batch_key=batch_key,
                    counts_layer=counts_layer,
                    n_pcs=n_pcs,
                    n_neighbors=n_neighbors,
                    random_state=random_state,
                )
            )
    return pd.DataFrame(rows)


def inverse_average_rank(
    results: pd.DataFrame,
    *,
    metrics: Iterable[str] = ("ARI", "NMI", "kNN_acc", "ASW", "VarRatio", "invLISI"),
) -> pd.DataFrame:
    """Rank methods within each dataset/K and return inverse average rank."""
    frames = []
    metrics = [metric for metric in metrics if metric in results.columns]
    for (dataset, k), sub in results.groupby(["dataset", "K"], dropna=False):
        ranked = sub[["dataset", "method", "K"]].copy()
        for metric in metrics:
            ranked[f"{metric}_rank"] = sub[metric].rank(ascending=False, method="average")
        rank_cols = [f"{metric}_rank" for metric in metrics]
        ranked["average_rank"] = ranked[rank_cols].mean(axis=1)
        ranked["inverse_average_rank"] = 1.0 / ranked["average_rank"]
        frames.append(ranked)
    return pd.concat(frames, ignore_index=True)
