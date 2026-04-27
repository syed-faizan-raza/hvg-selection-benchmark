from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from .metrics import inverse_average_rank


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
