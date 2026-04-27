from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .hvg import METHOD_ORDER, select_top_genes
from .metrics import inverse_average_rank, prepare_embedding


METHOD_COLORS = {
    "Seurat v3": "#4E79A7",
    "Seurat": "#F28E2B",
    "scran MV": "#E15759",
    "SCT proxy": "#76B7B2",
}


def save_metric_lines(
    results: pd.DataFrame,
    *,
    metric: str,
    output: str | Path,
    title: str | None = None,
):
    """Line plot of one metric across K, faceted by dataset."""
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)

    g = sns.relplot(
        data=results,
        x="K",
        y=metric,
        hue="method",
        col="dataset",
        kind="line",
        marker="o",
        palette=METHOD_COLORS,
        height=3.2,
        aspect=1.25,
        facet_kws={"sharey": False},
    )
    g.set_axis_labels("# HVGs (K)", metric)
    if title:
        g.figure.suptitle(title, y=1.05)
    g.figure.savefig(output, dpi=350, bbox_inches="tight")
    plt.close(g.figure)


def save_rank_summary(results: pd.DataFrame, *, output: str | Path):
    """Plot inverse average rank across all configured metrics."""
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    ranks = inverse_average_rank(results)

    fig, ax = plt.subplots(figsize=(8, 3.8))
    sns.lineplot(
        data=ranks,
        x="K",
        y="inverse_average_rank",
        hue="method",
        style="dataset",
        markers=True,
        dashes=False,
        palette=METHOD_COLORS,
        ax=ax,
    )
    ax.set_xlabel("# HVGs (K)")
    ax.set_ylabel("Inverse average rank")
    ax.set_title("Performance ranking across clustering metrics")
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(output, dpi=350, bbox_inches="tight")
    plt.close(fig)
    return ranks


def save_umap_grid(
    adata,
    scores_by_method: dict[str, pd.Series],
    *,
    k_values: Iterable[int],
    output: str | Path,
    label_key: str | None = None,
    counts_layer: str = "counts",
    n_pcs: int = 30,
    n_neighbors: int = 15,
    random_state: int = 7,
    method_order: Iterable[str] = METHOD_ORDER,
):
    """Create a method-by-K UMAP grid from selected HVGs."""
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)

    methods = [method for method in method_order if method in scores_by_method]
    k_values = [int(k) for k in k_values]
    fig, axes = plt.subplots(
        len(methods),
        len(k_values),
        figsize=(2.8 * len(k_values), 2.55 * len(methods)),
        squeeze=False,
    )

    labels = None
    palette = None
    if label_key and label_key in adata.obs:
        labels = adata.obs[label_key].astype(str)
        categories = sorted(labels.unique())
        colors = sns.color_palette("tab20", n_colors=len(categories))
        palette = dict(zip(categories, colors, strict=False))

    for row, method in enumerate(methods):
        for col, k in enumerate(k_values):
            ax = axes[row, col]
            genes = select_top_genes(scores_by_method[method], k)
            ad_emb = prepare_embedding(
                adata,
                genes,
                counts_layer=counts_layer,
                n_pcs=n_pcs,
                n_neighbors=n_neighbors,
                random_state=random_state,
                compute_umap=True,
            )
            coords = ad_emb.obsm["X_umap"]

            if labels is not None and palette is not None:
                colors = labels.loc[ad_emb.obs_names].map(palette)
                ax.scatter(coords[:, 0], coords[:, 1], s=4, c=list(colors), linewidths=0, alpha=0.85)
            else:
                ax.scatter(coords[:, 0], coords[:, 1], s=4, color="#4E79A7", linewidths=0, alpha=0.85)

            if row == 0:
                ax.set_title(f"K={k}", fontsize=10)
            if col == 0:
                ax.set_ylabel(method, fontsize=10)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.spines[["top", "right", "bottom", "left"]].set_visible(False)

    if labels is not None and palette is not None:
        handles = [
            plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=color, markersize=5, label=cat)
            for cat, color in palette.items()
        ]
        fig.legend(
            handles=handles,
            loc="lower center",
            ncol=min(6, len(handles)),
            frameon=False,
            fontsize=8,
            bbox_to_anchor=(0.5, -0.01),
        )
        bottom = 0.08
    else:
        bottom = 0.02

    fig.tight_layout(rect=[0, bottom, 1, 1])
    fig.savefig(output, dpi=350, bbox_inches="tight")
    plt.close(fig)

