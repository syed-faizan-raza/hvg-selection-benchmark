from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from statsmodels.nonparametric.smoothers_lowess import lowess

METHOD_ORDER = ("Seurat v3", "Seurat", "scran MV", "SCT proxy")


def _counts_matrix(adata, counts_layer: str = "counts"):
    if counts_layer in adata.layers:
        return adata.layers[counts_layer]
    return adata.X


def _as_clean_series(values, adata, name: str) -> pd.Series:
    values = np.asarray(values, dtype=float).ravel()
    finite = np.isfinite(values)
    high = float(np.nanmax(values[finite])) if finite.any() else 0.0
    values = np.nan_to_num(values, nan=0.0, posinf=high, neginf=0.0)
    return pd.Series(values, index=adata.var_names.astype(str), name=name)


def _pick_var_score(var, preferred: Iterable[str]) -> np.ndarray:
    for key in preferred:
        if key in var.columns:
            return var[key].fillna(0.0).to_numpy()
    raise KeyError(f"No expected HVG score column found. Tried: {', '.join(preferred)}")


def _set_x_to_counts(adata, counts_layer: str):
    if counts_layer in adata.layers:
        adata.X = adata.layers[counts_layer].copy()
    return adata


def _mean_var(X):
    if sparse.issparse(X):
        mean = np.asarray(X.mean(axis=0)).ravel()
        mean_sq = np.asarray(X.power(2).mean(axis=0)).ravel()
        var = np.maximum(mean_sq - mean**2, 0.0)
    else:
        X = np.asarray(X)
        mean = X.mean(axis=0)
        var = X.var(axis=0)
    return mean, var


def score_seurat_v3(adata, *, counts_layer: str = "counts") -> pd.Series:
    """Seurat v3 flavor on raw counts."""
    ad = adata.copy()
    _set_x_to_counts(ad, counts_layer)
    sc.pp.highly_variable_genes(
        ad,
        flavor="seurat_v3",
        layer=counts_layer if counts_layer in ad.layers else None,
        n_top_genes=ad.n_vars,
        subset=False,
        check_values=True,
    )
    score = _pick_var_score(ad.var, ("variances_norm", "variances"))
    return _as_clean_series(score, ad, "Seurat v3")


def score_seurat(adata, *, counts_layer: str = "counts") -> pd.Series:
    """Classic Seurat dispersion score after library-size normalization and log1p."""
    ad = adata.copy()
    _set_x_to_counts(ad, counts_layer)
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    sc.pp.highly_variable_genes(ad, flavor="seurat", n_top_genes=ad.n_vars, subset=False)
    score = _pick_var_score(ad.var, ("dispersions_norm", "dispersions", "variances_norm"))
    return _as_clean_series(score, ad, "Seurat")


def score_scran_mv(adata, *, counts_layer: str = "counts", lowess_frac: float = 0.25) -> pd.Series:
    """scran-style mean-variance residuals on log-normalized expression."""
    ad = adata.copy()
    _set_x_to_counts(ad, counts_layer)
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)

    mean, var = _mean_var(ad.X)
    x = np.log1p(mean)
    y = np.log1p(var)
    ok = np.isfinite(x) & np.isfinite(y) & (mean > 0)

    residual = np.zeros(ad.n_vars, dtype=float)
    if ok.sum() >= 5:
        order = np.argsort(x[ok])
        trend_sorted = lowess(
            y[ok][order],
            x[ok][order],
            frac=lowess_frac,
            return_sorted=False,
        )
        trend = np.interp(x, x[ok][order], trend_sorted)
        residual = y - trend

    return _as_clean_series(residual, ad, "scran MV")


def score_sct_proxy(
    adata,
    *,
    counts_layer: str = "counts",
    clip: float = 30.0,
    block_size: int = 2000,
) -> pd.Series:
    """SCT-like proxy: variance of clipped Pearson residuals from raw counts."""
    X = _counts_matrix(adata, counts_layer)
    n_cells, n_genes = X.shape

    library_size = np.asarray(X.sum(axis=1)).ravel().astype(float) + 1e-8
    gene_total = np.asarray(X.sum(axis=0)).ravel().astype(float)
    gene_prob = gene_total / (gene_total.sum() + 1e-12)

    scores = np.zeros(n_genes, dtype=float)
    for start in range(0, n_genes, block_size):
        stop = min(start + block_size, n_genes)
        X_block = X[:, start:stop]
        X_block = X_block.toarray() if sparse.issparse(X_block) else np.asarray(X_block)
        mu = library_size[:, None] * gene_prob[start:stop][None, :]
        resid = (X_block - mu) / np.sqrt(mu + 1.0)
        resid = np.clip(resid, -clip, clip)
        scores[start:stop] = np.mean(resid**2, axis=0)

    return _as_clean_series(scores, adata, "SCT proxy")


SCORERS = {
    "Seurat v3": score_seurat_v3,
    "Seurat": score_seurat,
    "scran MV": score_scran_mv,
    "SCT proxy": score_sct_proxy,
}


def score_all_methods(adata, methods: Iterable[str] | None = None, *, counts_layer: str = "counts"):
    """Return one per-gene score Series for each requested method."""
    methods = tuple(methods or METHOD_ORDER)
    scores = {}
    for method in methods:
        try:
            scorer = SCORERS[method]
        except KeyError as exc:
            available = ", ".join(SCORERS)
            raise KeyError(f"Unknown method '{method}'. Available: {available}") from exc
        scores[method] = scorer(adata, counts_layer=counts_layer)
    return scores


def select_top_genes(scores: pd.Series, k: int) -> list[str]:
    """Select the top K genes from a score vector."""
    k = min(int(k), len(scores))
    return scores.sort_values(ascending=False).index[:k].astype(str).tolist()


def scores_to_long_table(scores_by_method: dict[str, pd.Series], *, dataset: str) -> pd.DataFrame:
    """Convert method score Series into a tidy table with ranks."""
    frames = []
    for method, scores in scores_by_method.items():
        table = (
            scores.sort_values(ascending=False)
            .rename("score")
            .reset_index()
            .rename(columns={"index": "gene"})
        )
        table.insert(0, "dataset", dataset)
        table.insert(1, "method", method)
        table["rank"] = np.arange(1, len(table) + 1)
        frames.append(table)
    return pd.concat(frames, ignore_index=True)

