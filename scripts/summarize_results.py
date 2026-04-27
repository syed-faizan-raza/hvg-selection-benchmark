#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from hvg_benchmark.metrics import inverse_average_rank
from hvg_benchmark.plotting import save_metric_lines, save_rank_summary


def parse_args():
    parser = argparse.ArgumentParser(description="Summarize clustering benchmark results.")
    parser.add_argument("--metrics", required=True, help="CSV from run_clustering_benchmark.py")
    parser.add_argument("--outdir", default="outputs/summary")
    return parser.parse_args()


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    results = pd.read_csv(args.metrics)
    ranks = inverse_average_rank(results)
    ranks_path = outdir / "inverse_average_ranks.csv"
    ranks.to_csv(ranks_path, index=False)

    for metric in ["ARI", "NMI", "kNN_acc", "ASW", "VarRatio", "invLISI"]:
        if metric in results.columns:
            save_metric_lines(results, metric=metric, output=outdir / f"{metric}_by_K.png")
    save_rank_summary(results, output=outdir / "inverse_average_rank_by_K.png")

    print(f"Saved rank table: {ranks_path}")
    print(f"Saved summary figures under: {outdir}")


if __name__ == "__main__":
    main()
