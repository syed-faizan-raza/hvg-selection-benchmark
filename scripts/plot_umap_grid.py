#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from hvg_benchmark.hvg import METHOD_ORDER, score_all_methods
from hvg_benchmark.io import load_configured_dataset
from hvg_benchmark.plotting import save_umap_grid


def parse_args():
    parser = argparse.ArgumentParser(description="Plot stacked UMAPs for method x K comparisons.")
    parser.add_argument("--config", default="configs/datasets.local.yaml")
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--scores", default=None, help="Optional CSV from run_hvg_scores.py")
    parser.add_argument("--output", default=None)
    parser.add_argument("--k", nargs="+", type=int, default=None)
    parser.add_argument("--methods", nargs="+", default=None)
    return parser.parse_args()


def load_scores(path: str | Path):
    table = pd.read_csv(path)
    return {
        method: sub.set_index("gene")["score"].sort_values(ascending=False)
        for method, sub in table.groupby("method")
    }


def main():
    args = parse_args()
    adata, ds, config = load_configured_dataset(args.config, args.dataset)
    settings = config.get("benchmark", {})

    k_values = args.k or settings.get("k_values", [500, 1000, 2000])
    methods = args.methods or settings.get("methods", list(METHOD_ORDER))
    n_pcs = int(settings.get("n_pcs", 30))
    n_neighbors = int(settings.get("n_neighbors", 15))
    random_state = int(settings.get("random_state", 7))

    scores = load_scores(args.scores) if args.scores else score_all_methods(adata, methods=methods)
    output = Path(args.output or f"outputs/{ds.name}_umap_method_by_k.png")

    save_umap_grid(
        adata,
        scores,
        k_values=k_values,
        output=output,
        label_key=ds.label_key,
        n_pcs=n_pcs,
        n_neighbors=n_neighbors,
        random_state=random_state,
        method_order=methods,
    )
    print(f"Saved UMAP grid: {output}")


if __name__ == "__main__":
    main()
