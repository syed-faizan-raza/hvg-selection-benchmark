#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from hvg_benchmark.hvg import METHOD_ORDER, score_all_methods
from hvg_benchmark.io import load_configured_dataset
from hvg_benchmark.metrics import evaluate_scores_grid


def parse_args():
    parser = argparse.ArgumentParser(description="Run the method x K clustering benchmark.")
    parser.add_argument("--config", default="configs/datasets.local.yaml")
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--output", default=None)
    parser.add_argument("--k", nargs="+", type=int, default=None)
    parser.add_argument("--methods", nargs="+", default=None)
    return parser.parse_args()


def main():
    args = parse_args()
    adata, ds, config = load_configured_dataset(args.config, args.dataset)
    settings = config.get("benchmark", {})

    k_values = args.k or settings.get("k_values", [200, 500, 1000, 1500, 2000, 3000])
    methods = args.methods or settings.get("methods", list(METHOD_ORDER))
    n_pcs = int(settings.get("n_pcs", 30))
    n_neighbors = int(settings.get("n_neighbors", 15))
    random_state = int(settings.get("random_state", 7))

    if ds.label_key not in adata.obs:
        raise KeyError(
            f"Configured label_key '{ds.label_key}' is not present in adata.obs. "
            "Add labels or edit configs/datasets.local.yaml."
        )

    print(f"Loaded {ds.name}: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    print(f"Methods: {', '.join(methods)}")
    print(f"K values: {k_values}")

    scores = score_all_methods(adata, methods=methods)
    results = evaluate_scores_grid(
        adata,
        scores,
        dataset=ds.name,
        k_values=k_values,
        label_key=ds.label_key,
        batch_key=ds.batch_key,
        n_pcs=n_pcs,
        n_neighbors=n_neighbors,
        random_state=random_state,
    )

    output = Path(args.output or f"outputs/{ds.name}_clustering_metrics.csv")
    output.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(output, index=False)
    print(f"Saved clustering metrics: {output}")


if __name__ == "__main__":
    main()
