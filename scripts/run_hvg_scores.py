#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from hvg_benchmark.hvg import METHOD_ORDER, score_all_methods, scores_to_long_table
from hvg_benchmark.io import load_configured_dataset


def parse_args():
    parser = argparse.ArgumentParser(description="Score genes with each HVG method.")
    parser.add_argument("--config", default="configs/datasets.local.yaml")
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--output", default=None)
    parser.add_argument("--top-n", type=int, default=100)
    parser.add_argument("--methods", nargs="+", default=list(METHOD_ORDER))
    return parser.parse_args()


def main():
    args = parse_args()
    adata, ds, _ = load_configured_dataset(args.config, args.dataset)
    print(f"Loaded {ds.name}: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    scores = score_all_methods(adata, methods=args.methods)
    table = scores_to_long_table(scores, dataset=ds.name)

    output = Path(args.output or f"outputs/{ds.name}_hvg_scores.csv")
    output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output, index=False)

    top_output = output.with_name(output.stem.replace("_hvg_scores", "_top_genes") + ".csv")
    table.query("rank <= @args.top_n").to_csv(top_output, index=False)

    print(f"Saved full scores: {output}")
    print(f"Saved top genes:   {top_output}")


if __name__ == "__main__":
    main()
