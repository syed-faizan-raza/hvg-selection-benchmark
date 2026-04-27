#!/usr/bin/env python3
from __future__ import annotations

import argparse
import itertools
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


METRICS = ["ARI", "NMI", "kNN_acc", "ASW", "VarRatio", "invLISI"]
METHOD_ORDER = ["Seurat v3", "Seurat", "scran MV", "SCT proxy"]
METHOD_COLORS = {
    "Seurat v3": "#4E79A7",
    "Seurat": "#F28E2B",
    "scran MV": "#E15759",
    "SCT proxy": "#76B7B2",
}


def parse_args():
    parser = argparse.ArgumentParser(description="Build compact example report tables and figures.")
    parser.add_argument(
        "--metrics",
        default="docs/example_results/full_dataset_clustering_metrics.csv",
        help="Clustering metrics CSV produced by run_clustering_benchmark.py",
    )
    parser.add_argument(
        "--top-genes",
        default="docs/example_results/full_dataset_top_genes.csv",
        help="Top-gene CSV produced by run_hvg_scores.py",
    )
    parser.add_argument("--outdir", default="docs/example_results")
    return parser.parse_args()


def write_rank_tables(metrics: pd.DataFrame, outdir: Path):
    rank_rows = []
    for (dataset, k), sub in metrics.groupby(["dataset", "K"]):
        sub = sub.copy()
        for metric in METRICS:
            sub[f"{metric}_rank"] = sub[metric].rank(ascending=False, method="average")
        rank_cols = [f"{metric}_rank" for metric in METRICS]
        sub["average_rank"] = sub[rank_cols].mean(axis=1)
        sub["inverse_average_rank"] = 1 / sub["average_rank"]
        rank_rows.append(sub[["dataset", "method", "K", "average_rank", "inverse_average_rank"]])
    ranks = pd.concat(rank_rows, ignore_index=True)
    ranks.to_csv(outdir / "full_dataset_rank_summary.csv", index=False)

    best_rows = []
    for (dataset, k), sub in metrics.groupby(["dataset", "K"]):
        for metric in METRICS:
            idx = sub[metric].idxmax()
            best_rows.append(
                {
                    "dataset": dataset,
                    "K": k,
                    "metric": metric,
                    "best_method": sub.loc[idx, "method"],
                    "best_value": sub.loc[idx, metric],
                }
            )
    pd.DataFrame(best_rows).to_csv(outdir / "full_dataset_best_methods.csv", index=False)
    return ranks


def write_overlap_table(top_genes: pd.DataFrame, outdir: Path, top_n: int = 100):
    methods = [method for method in METHOD_ORDER if method in set(top_genes["method"])]
    gene_sets = {
        method: set(top_genes[(top_genes.method == method) & (top_genes["rank"] <= top_n)]["gene"])
        for method in methods
    }

    rows = []
    for method_a, method_b in itertools.combinations(methods, 2):
        shared = gene_sets[method_a] & gene_sets[method_b]
        union = gene_sets[method_a] | gene_sets[method_b]
        rows.append(
            {
                "method_a": method_a,
                "method_b": method_b,
                "top_n": top_n,
                "overlap_count": len(shared),
                "jaccard": len(shared) / len(union) if union else np.nan,
                "shared_genes": ";".join(sorted(shared)),
            }
        )
    pd.DataFrame(rows).to_csv(outdir / "full_dataset_top_gene_overlap.csv", index=False)


def save_top_gene_panel(top_genes: pd.DataFrame, outdir: Path):
    methods = [method for method in METHOD_ORDER if method in set(top_genes["method"])]
    fig, axes = plt.subplots(1, len(methods), figsize=(13.5, 4.2), sharex=False)

    for ax, method in zip(axes, methods):
        sub = top_genes[top_genes.method == method].sort_values("rank").head(10).iloc[::-1]
        scaled = sub["score"].astype(float)
        scaled = scaled / scaled.max() if scaled.max() else scaled
        ax.barh(sub["gene"], scaled, color=METHOD_COLORS[method], alpha=0.9)
        ax.set_title(method, fontsize=11, pad=8)
        ax.set_xlim(0, 1.08)
        ax.set_xlabel("scaled score", fontsize=9)
        ax.tick_params(axis="y", labelsize=8)
        ax.tick_params(axis="x", labelsize=8)
        ax.grid(axis="x", alpha=0.22)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.suptitle("Top HVGs by method, full-dataset example run", fontsize=13, y=1.03)
    fig.subplots_adjust(left=0.08, right=0.99, top=0.83, bottom=0.17, wspace=0.55)
    fig.savefig(outdir / "full_dataset_top_genes_panel.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_rank_plot(ranks: pd.DataFrame, outdir: Path):
    fig, ax = plt.subplots(figsize=(7.5, 4.3))
    for method in [method for method in METHOD_ORDER if method in set(ranks["method"])]:
        sub = ranks[ranks.method == method].sort_values("K")
        ax.plot(
            sub["K"],
            sub["inverse_average_rank"],
            marker="o",
            linewidth=2,
            color=METHOD_COLORS[method],
            label=method,
        )

    ax.set_xticks(sorted(ranks["K"].unique()))
    ax.set_xlabel("HVG count (K)")
    ax.set_ylabel("inverse average rank")
    ax.set_title("Average performance rank across metrics")
    ax.grid(axis="y", alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(outdir / "full_dataset_rank_summary.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    metrics = pd.read_csv(args.metrics)
    top_genes = pd.read_csv(args.top_genes)
    ranks = write_rank_tables(metrics, outdir)
    write_overlap_table(top_genes, outdir)
    save_top_gene_panel(top_genes, outdir)
    save_rank_plot(ranks, outdir)
    print(f"Saved example report files under: {outdir}")


if __name__ == "__main__":
    main()

