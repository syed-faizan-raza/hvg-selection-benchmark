# HVG Selection Benchmark

Companion code for the poster:

**Systematic Evaluation of Highly Variable Gene Selection and Batch-Aware Approaches in Single-Cell RNA-seq**

Authors: Syed Faizan, Binte Zehra, Bakhrom K. Berdiev

This repository is intentionally lightweight. It contains the poster PDF and the scripts used to organize the HVG-selection benchmark workflow. Large datasets and generated result tables are not included.

## Poster

[docs/poster/SyedFaizan_Poster.pdf](docs/poster/SyedFaizan_Poster.pdf)

## What This Code Covers

The scripts compare four highly variable gene selection methods:

- Seurat
- Seurat v3
- scran-MV
- SCT proxy

Each method can be evaluated across different numbers of selected genes (`K`) using the same downstream workflow:

```text
HVG selection -> normalization/log1p -> PCA -> neighbors -> Leiden clustering
```

The intended clustering metrics are:

```text
ARI, NMI, kNN accuracy, ASW, variance ratio, invLISI
```

If a batch column is available, the scripts can also report batch LISI.

## Folder Structure

```text
.
├── configs/          # Example dataset configuration
├── data/             # Local data location, ignored by git
├── docs/poster/      # Poster PDF
├── scripts/          # Run scripts
└── src/              # Shared helper code used by the scripts
```

## Setup

Install the usual Scanpy stack in your preferred environment. The scripts expect:

```text
scanpy
anndata
numpy
pandas
scipy
scikit-learn
statsmodels
matplotlib
seaborn
pyyaml
leidenalg
igraph
scikit-misc
```

Or install the listed dependencies with:

```bash
pip install -r requirements.txt
```

## Data

Raw data are not committed.

Copy the example config:

```bash
cp configs/datasets.example.yaml configs/datasets.local.yaml
```

Then edit `configs/datasets.local.yaml` so each dataset points to a local `.h5ad`, 10x `.h5`, or 10x matrix directory.

For supervised clustering metrics, the `AnnData.obs` table must contain the configured `label_key`.

## Scripts

Score genes:

```bash
python scripts/run_hvg_scores.py \
  --config configs/datasets.local.yaml \
  --dataset pbmc10k
```

Run clustering metrics:

```bash
python scripts/run_clustering_benchmark.py \
  --config configs/datasets.local.yaml \
  --dataset pbmc10k \
  --k 200 500 1000 1500 2000 3000 6000
```

Plot a method-by-K UMAP grid:

```bash
python scripts/plot_umap_grid.py \
  --config configs/datasets.local.yaml \
  --dataset pbmc10k \
  --k 500 1000 2000
```

Generated files go to `outputs/` by default, which is ignored by git.

## Example Output

A full-dataset documentation output folder is included so the expected file shapes are visible without rerunning the analysis:

- [docs/example_results/full_dataset_clustering_metrics.csv](docs/example_results/full_dataset_clustering_metrics.csv)
- [docs/example_results/full_dataset_umap_grid.png](docs/example_results/full_dataset_umap_grid.png)

These outputs were generated from the full labelled dataset available locally in this project workspace (`10,942` cells after loading; `12,205` genes after the 1% detected-cell filter) using `K = 500, 1000, 2000`. They are included as a real documentation run, not as the final two-dataset PBMC/CHD poster rerun.

## Current Status

This is a cleaned-up companion repository for the poster, not a fully rerun archive. The included full-dataset documentation run shows the expected workflow output. The remaining final piece is to rerun the complete PBMC/CHD benchmark locally and generate one combined table:

```text
dataset, method, K, ARI, NMI, kNN_acc, ASW, VarRatio, invLISI
```

Once that table is regenerated, the poster's PBMC-versus-CHD clustering statements can be reproduced directly from code.
