#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


METRICS = [
    ("ARI", "ARI"),
    ("NMI", "NMI"),
    ("kNN_acc", "kNN accuracy"),
    ("ASW", "ASW"),
    ("VarRatio", "Variance ratio"),
    ("invLISI", "invLISI"),
]

METHOD_COLORS = {
    "Seurat v3": "#4E79A7",
    "Seurat": "#F28E2B",
    "scran MV": "#E15759",
    "SCT proxy": "#76B7B2",
}


def parse_args():
    parser = argparse.ArgumentParser(description="Plot benchmark metric panels from a CSV file.")
    parser.add_argument(
        "--metrics",
        default="docs/example_results/full_dataset_clustering_metrics.csv",
        help="CSV produced by run_clustering_benchmark.py",
    )
    parser.add_argument(
        "--output",
        default="docs/example_results/full_dataset_metric_panels.png",
        help="Output image path",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    df = pd.read_csv(args.metrics)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    methods = [m for m in METHOD_COLORS if m in set(df["method"])]
    remaining = [m for m in df["method"].drop_duplicates() if m not in methods]
    methods.extend(remaining)
    k_values = sorted(df["K"].unique())

    fig, axes = plt.subplots(2, 3, figsize=(11.2, 6.8), sharex=True)
    axes = axes.ravel()

    for ax, (metric, label) in zip(axes, METRICS):
        for method in methods:
            sub = df[df["method"] == method].sort_values("K")
            if sub.empty or metric not in sub:
                continue
            ax.plot(
                sub["K"],
                sub[metric],
                marker="o",
                linewidth=2.0,
                markersize=4.8,
                color=METHOD_COLORS.get(method),
                label=method,
            )
        ax.set_title(label, fontsize=11, pad=9)
        ax.set_xlabel("HVG count (K)", fontsize=9)
        ax.set_xticks(k_values)
        ax.grid(axis="y", alpha=0.25, linewidth=0.8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(axis="both", labelsize=8)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        ncol=min(4, len(labels)),
        frameon=False,
        fontsize=10,
        bbox_to_anchor=(0.5, 0.955),
    )
    fig.suptitle("Full-dataset clustering metrics by HVG method", fontsize=13, y=0.995)
    fig.subplots_adjust(left=0.07, right=0.98, top=0.80, bottom=0.08, hspace=0.52, wspace=0.35)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved metric panels: {output}")


if __name__ == "__main__":
    main()
