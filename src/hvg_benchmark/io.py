from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import scanpy as sc
import yaml


@dataclass(frozen=True)
class DatasetConfig:
    name: str
    source: str
    label_key: str | None = None
    batch_key: str | None = None
    description: str = ""
    min_cells_fraction: float = 0.0


def read_config(path: str | Path) -> dict[str, Any]:
    """Read a YAML benchmark configuration."""
    with Path(path).open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def dataset_config(config: dict[str, Any], name: str) -> DatasetConfig:
    """Return a normalized dataset config from the YAML dictionary."""
    try:
        raw = config["datasets"][name]
    except KeyError as exc:
        available = ", ".join(sorted(config.get("datasets", {})))
        raise KeyError(f"Dataset '{name}' not found. Available: {available}") from exc

    return DatasetConfig(
        name=name,
        source=str(raw["source"]),
        label_key=raw.get("label_key"),
        batch_key=raw.get("batch_key"),
        description=raw.get("description", ""),
        min_cells_fraction=float(raw.get("min_cells_fraction", 0.0)),
    )


def ensure_counts_layer(adata, layer: str = "counts"):
    """Keep a raw-count layer for methods that must run before normalization."""
    if layer not in adata.layers:
        adata.layers[layer] = adata.X.copy()
    return adata


def load_adata(source: str | Path, *, counts_layer: str = "counts"):
    """Load an AnnData object from a Scanpy-supported source."""
    source = str(source)
    lower = source.lower()

    if source == "pbmc3k":
        adata = sc.datasets.pbmc3k()
    elif lower.endswith(".h5ad"):
        adata = sc.read_h5ad(source)
    elif lower.endswith(".h5"):
        adata = sc.read_10x_h5(source)
    else:
        path = Path(source)
        if path.is_dir():
            adata = sc.read_10x_mtx(path)
        else:
            raise ValueError(f"Unsupported data source: {source}")

    adata.var_names_make_unique()
    return ensure_counts_layer(adata, counts_layer)


def filter_counts(adata, *, min_cells_fraction: float = 0.0, counts_layer: str = "counts"):
    """Drop empty cells and genes detected in fewer than the requested cells."""
    import scipy.sparse as sp

    X = adata.layers[counts_layer] if counts_layer in adata.layers else adata.X

    cell_counts = np.asarray(X.sum(axis=1)).ravel()
    if np.any(cell_counts <= 0):
        adata = adata[cell_counts > 0].copy()
        X = adata.layers[counts_layer] if counts_layer in adata.layers else adata.X

    if min_cells_fraction > 0:
        detected = np.asarray((X > 0).sum(axis=0)).ravel() if sp.issparse(X) else (X > 0).sum(axis=0)
        min_cells = max(1, int(round(min_cells_fraction * adata.n_obs)))
        adata = adata[:, detected >= min_cells].copy()

    return ensure_counts_layer(adata, counts_layer)


def load_configured_dataset(config_path: str | Path, dataset: str, *, counts_layer: str = "counts"):
    """Load and lightly filter a dataset described by the YAML config."""
    config = read_config(config_path)
    ds = dataset_config(config, dataset)
    adata = load_adata(ds.source, counts_layer=counts_layer)
    adata = filter_counts(
        adata,
        min_cells_fraction=ds.min_cells_fraction,
        counts_layer=counts_layer,
    )
    return adata, ds, config
