"""Marker-centric plots mirroring the R featureplot & heatmap workflow."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, List

import anndata as ad
import pandas as pd
import scanpy as sc

from . import config
from .utils.io import ensure_dirs
from .utils.logging import get_logger
from .utils.markers import BASIC_MARKER_PANEL

LOGGER = get_logger(__name__)


FEATURE_SETS = {
    "neutrophil": ["GOS2", "S100A9", "S100A8", "CXCL8"],
    "plasma": ["IGHG1", "MZB1", "SDC1", "CD79A"],
}


def load_celltype_data(path: Path | None = None) -> ad.AnnData:
    path = path or (config.DEFAULT_OUTPUT_SUBDIRS["celltype"] / "sce_celltype.h5ad")
    if not path.exists():
        raise FileNotFoundError(f"Annotated AnnData not found: {path}")
    adata = ad.read_h5ad(path)
    # Ensure obs_names are unique to avoid indexing issues
    if not adata.obs_names.is_unique:
        adata.obs_names_make_unique()
    return adata


def highlight_clusters(adata: ad.AnnData, clusters: Iterable[int], *, label: str) -> None:
    mask = adata.obs["seurat_clusters"].isin(clusters)
    color = pd.Series("grey", index=adata.obs_names)
    color[mask] = "red"
    adata.obs[f"highlight_{label}"] = color
    sc.pl.tsne(
        adata,
        color=f"highlight_{label}",
        palette=["grey", "red"],
        size=20,
        save=f"_{label}_highlight.png",
        show=False,
    )


def feature_plot_panel(adata: ad.AnnData, genes: List[str], *, tag: str) -> None:
    # Filter genes to only include those present in the dataset
    available_genes = [gene for gene in genes if gene in adata.var_names]
    if available_genes:
        sc.pl.tsne(
            adata,
            color=available_genes,
            cmap="Reds",
            ncols=2,
            save=f"_{tag}_features.png",
            show=False,
        )
    else:
        LOGGER.warning(f"No genes from {tag} panel found in dataset")


def blended_plot(adata: ad.AnnData, genes: List[str], *, tag: str) -> None:
    if len(genes) != 2:
        LOGGER.warning("Blend plot expects exactly two genes; received %s", genes)
        return
    # Filter genes to only include those present in the dataset
    available_genes = [gene for gene in genes if gene in adata.var_names]
    if len(available_genes) == 2:
        sc.pl.tsne(
            adata,
            color=available_genes,
            color_map="viridis",
            blend=True,
            save=f"_{tag}_blend.png",
            show=False,
        )
    else:
        LOGGER.warning(f"Not enough genes from {tag} blend found in dataset: {available_genes}")


def compute_top_markers(adata: ad.AnnData, n_top: int = 5) -> pd.DataFrame:
    sc.tl.rank_genes_groups(
        adata,
        groupby="celltype",
        method="wilcoxon",
        n_genes=n_top,
        use_raw=False,
    )
    return sc.get.rank_genes_groups_df(adata, None)


def plot_marker_heatmap(adata: ad.AnnData, marker_df: pd.DataFrame, *, n_cells: int = 100) -> None:
    # Filter out rows with NaN values that might cause issues
    marker_df = marker_df.dropna(subset=["names", "group"])

    if marker_df.empty:
        LOGGER.warning("No valid markers found for heatmap")
        return

    marker_map = (
        marker_df.groupby("group")["names"].apply(list).to_dict()
    )
    genes = sorted({gene for genes in marker_map.values() for gene in genes})

    # Filter genes to only those present in the dataset
    genes = [gene for gene in genes if gene in adata.var_names]

    if not genes:
        LOGGER.warning("No marker genes found in dataset for heatmap")
        return

    # Sample cells from each celltype, handling non-unique indices
    sampled_indices = []
    for celltype, group in adata.obs.groupby("celltype"):
        sampled_indices.extend(group.head(n_cells).index.tolist())
    subset = adata[sampled_indices, :]
    sc.pl.heatmap(
        subset,
        var_names=genes,
        groupby="celltype",
        swap_axes=True,
        standard_scale="var",
        save="_markers_heatmap.png",
        show=False,
    )


def plot_marker_dotplot(adata: ad.AnnData, marker_df: pd.DataFrame) -> None:
    genes = marker_df["names"].dropna().unique().tolist()
    # Filter genes to only those present in the dataset
    genes = [gene for gene in genes if gene in adata.var_names]

    if not genes:
        LOGGER.warning("No marker genes found in dataset for dotplot")
        return

    sc.pl.dotplot(
        adata,
        var_names=genes,
        groupby="celltype",
        standard_scale="var",
        save="_markers_dotplot.png",
        show=False,
    )


def run_step4_feature_heatmap() -> None:
    out_dir = config.DEFAULT_OUTPUT_SUBDIRS["plots"]
    ensure_dirs([out_dir])
    sc.settings.figdir = str(out_dir)

    adata = load_celltype_data()

    highlight_clusters(adata, clusters=[4], label="neutrophil")
    feature_plot_panel(adata, FEATURE_SETS["neutrophil"], tag="neutrophil")
    highlight_clusters(adata, clusters=[7], label="plasma")
    feature_plot_panel(adata, FEATURE_SETS["plasma"], tag="plasma")

    feature_plot_panel(adata, BASIC_MARKER_PANEL, tag="marker_grid")
    # Note: blended_plot removed due to compatibility issues with current scanpy version
    # blended_plot(adata, ["S100A9", "S100A8"], tag="s100_blend")

    markers = compute_top_markers(adata)
    markers.to_csv(out_dir / "celltype_markers.csv", index=False)
    plot_marker_heatmap(adata, markers)
    plot_marker_dotplot(adata, markers)

