"""Cell clustering, marker detection, and annotation utilities."""

from __future__ import annotations

from pathlib import Path
from typing import Dict

import anndata as ad
import pandas as pd
import scanpy as sc

from . import config
from .utils.io import ensure_dirs
from .utils.logging import get_logger, time_block
from .utils.markers import BASIC_MARKER_PANEL, CELLTYPE_MAPPING

LOGGER = get_logger(__name__)


def load_harmony_object(path: Path | None = None) -> ad.AnnData:
    path = path or (config.DEFAULT_OUTPUT_SUBDIRS["harmony"] / "sce_all_int.h5ad")
    if not path.exists():
        raise FileNotFoundError(f"Integrated AnnData not found: {path}")
    LOGGER.info("Loading AnnData from %s", path)
    return ad.read_h5ad(path)


def pick_resolution(adata: ad.AnnData, resolution: float = 0.5) -> str:
    key = f"leiden_{resolution}"
    if key not in adata.obs:
        available = [col for col in adata.obs.columns if col.startswith("leiden_")]
        raise KeyError(f"Resolution {resolution} not found. Existing keys: {available}")
    adata.obs["seurat_clusters"] = adata.obs[key].astype("category")
    LOGGER.info("Selected Leiden clustering at resolution %.2f", resolution)
    return key


def annotate_celltypes(adata: ad.AnnData, mapping: Dict[int, str] | None = None) -> None:
    mapping = mapping or CELLTYPE_MAPPING
    # Convert seurat_clusters to int for mapping
    cluster_int = adata.obs["seurat_clusters"].astype(int)
    adata.obs["celltype"] = cluster_int.map(mapping).fillna("Unknown")


def compute_marker_statistics(adata: ad.AnnData) -> pd.DataFrame:
    LOGGER.info("Computing marker genes using rank_genes_groups")
    sc.tl.rank_genes_groups(
        adata,
        groupby="seurat_clusters",
        method="wilcoxon",
        n_genes=50,
        use_raw=False,
    )
    return sc.get.rank_genes_groups_df(adata, None)


def export_celltype_ann_data(adata: ad.AnnData, out_dir: Path | None = None) -> Path:
    out_dir = out_dir or config.DEFAULT_OUTPUT_SUBDIRS["celltype"]
    ensure_dirs([out_dir])
    output_path = out_dir / "sce_celltype.h5ad"
    adata.write_h5ad(output_path)
    LOGGER.info("Saved annotated AnnData to %s", output_path)
    return output_path


def run_step3_celltype(resolution: float = 0.5) -> None:
    out_dir = config.DEFAULT_OUTPUT_SUBDIRS["celltype"]
    ensure_dirs([out_dir])
    sc.settings.figdir = str(out_dir)

    timings_file = config.DEFAULT_OUTPUT_SUBDIRS["timings"] / "step03_celltype.jsonl"

    with time_block("load_harmony_object", write_jsonl=timings_file):
        adata = load_harmony_object()
    with time_block("pick_resolution", write_jsonl=timings_file):
        pick_resolution(adata, resolution=resolution)
    with time_block("annotate_celltypes", write_jsonl=timings_file):
        annotate_celltypes(adata)
    with time_block("compute_marker_statistics", write_jsonl=timings_file):
        markers = compute_marker_statistics(adata)
    markers_path = out_dir / "markers_rank_genes.csv"
    with time_block("write_markers_csv", write_jsonl=timings_file):
        markers.to_csv(markers_path, index=False)
    LOGGER.info("Exported marker table to %s", markers_path)
    with time_block("plot_umap", write_jsonl=timings_file):
        sc.pl.umap(
            adata,
            color=["seurat_clusters", "celltype", "orig.ident", "patient"],
            wspace=0.4,
            save="_step3_overview.png",
            show=False,
        )
    with time_block("plot_marker_dotplot", write_jsonl=timings_file):
        # Filter markers to only include genes present in the dataset
        available_markers = [gene for gene in BASIC_MARKER_PANEL if gene in adata.var_names]
        if available_markers:
            sc.pl.dotplot(
                adata,
                var_names=available_markers,
                groupby="celltype",
                standard_scale="var",
                dendrogram=False,
                save="_marker_panel.png",
                show=False,
            )
        else:
            LOGGER.warning("No marker genes found in dataset for dotplot")

    with time_block("export_celltype_ann_data", write_jsonl=timings_file):
        export_celltype_ann_data(adata, out_dir=out_dir)

