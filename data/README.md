# Data

Raw single-cell count matrices are not tracked in this repository.

Place the local datasets here, then copy `configs/datasets.example.yaml` to
`configs/datasets.local.yaml` and point each dataset entry to the correct file.

Supported inputs:

- `.h5ad`
- 10x Genomics `.h5`
- 10x Genomics matrix directory readable by `scanpy.read_10x_mtx`
- the built-in keyword `pbmc3k` for quick local examples

For clustering metrics such as ARI, NMI, kNN accuracy, ASW, variance ratio, and
invLISI, the loaded `AnnData.obs` table must contain the configured `label_key`.
