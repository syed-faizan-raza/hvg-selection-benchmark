# Full-Dataset Example Run

These files are compact outputs from a local full-dataset run:

- `10,942` cells
- `12,205` genes after filtering genes detected in at least 1% of cells
- HVG counts: `K = 500, 1000, 2000`
- Methods: Seurat v3, Seurat, scran-MV, SCT proxy

Included files:

- `full_dataset_clustering_metrics.csv`: clustering metrics by method and K
- `full_dataset_metric_panels.png`: metric panels plotted from the metrics CSV
- `full_dataset_top_genes.csv`: top 100 genes per method
- `full_dataset_top_genes_panel.png`: top 10 genes per method
- `full_dataset_top_gene_overlap.csv`: pairwise top-100 gene overlap
- `full_dataset_rank_summary.csv`: inverse average rank across metrics
- `full_dataset_rank_summary.png`: rank summary figure
- `full_dataset_best_methods.csv`: best method per metric and K
- `full_dataset_umap_grid.png`: full-data UMAP grid from the HVG selections

The raw `.h5ad` object is not included.

