"""T-cell subset annotation using marker-based scoring."""

from __future__ import annotations

from collections import OrderedDict
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

from . import config
from .utils.io import ensure_dirs
from .utils.logging import get_logger
from .utils.markers import T_CELL_SIGNATURES

LOGGER = get_logger(__name__)

def load_celltype_data() -> ad.AnnData:
    path = config.DEFAULT_OUTPUT_SUBDIRS["celltype"] / "sce_celltype.h5ad"
    if not path.exists():
        raise FileNotFoundError("Run step03_celltype.py first to create sce_celltype.h5ad")
    return ad.read_h5ad(path)

def subset_t_cells(adata: ad.AnnData) -> ad.AnnData:
    subset = adata[adata.obs["celltype"] == "T"].copy()
    LOGGER.info("T-cell subset size: %d", subset.n_obs)
    return subset

def preprocess(adata: ad.AnnData) -> None:
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    try:
        has_negative = float(adata.X.min()) < 0
    except Exception:
        has_negative = False
    already_logged = bool(adata.uns.get("log1p")) or has_negative

    if not already_logged:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    else:
        LOGGER.info("T-cell subset appears log-normalized; skipping normalize_total/log1p")

    hv_mask = None
    if has_negative and "highly_variable" in adata.var:
        mask = adata.var["highly_variable"].to_numpy().astype(bool)
        if mask.sum() == 0:
            LOGGER.warning("Existing highly_variable mask empty; recomputing for T-cell subset")
            sc.pp.highly_variable_genes(adata, n_top_genes=config.N_HVG, subset=False, flavor="pearson_residuals")
            hv_mask = adata.var["highly_variable"].to_numpy().astype(bool)
        else:
            LOGGER.info("Using existing highly_variable mask with %d genes for T cells", int(mask.sum()))
            hv_mask = mask
    else:
        sc.pp.highly_variable_genes(adata, n_top_genes=config.N_HVG, subset=False, flavor="seurat_v3")
        hv_mask = adata.var["highly_variable"].to_numpy().astype(bool)

    if hv_mask is None or hv_mask.sum() == 0:
        LOGGER.warning("No highly variable genes detected for T cells; using all genes")
        hv_mask = np.ones(adata.n_vars, dtype=bool)

    adata.var["highly_variable"] = hv_mask
    sc.pp.scale(adata)
    if np.isnan(adata.X).any():
        LOGGER.warning("NaNs detected after scaling T cells; replacing with zeros")
        adata.X = np.nan_to_num(adata.X, nan=0.0, posinf=0.0, neginf=0.0)

    sc.tl.pca(adata, n_comps=config.N_PCS, use_highly_variable=True)
    sc.pp.neighbors(adata, n_pcs=15)
    sc.tl.umap(adata)
    sc.tl.tsne(adata, random_state=config.SEED_TSNE)
    sc.tl.leiden(adata, resolution=0.1, key_added="tcell_clusters")

def score_signatures(adata: ad.AnnData) -> list[str]:
    score_cols: list[str] = []
    for label, genes in T_CELL_SIGNATURES.items():
        score_name = f"score_{label}"
        sc.tl.score_genes(adata, gene_list=genes, score_name=score_name)
        score_cols.append(score_name)
    scores = adata.obs[score_cols].to_numpy()
    best_idx = np.argmax(scores, axis=1)
    labels = list(T_CELL_SIGNATURES)
    adata.obs["tcell_label"] = [labels[i] for i in best_idx]
    return score_cols

def cluster_signature_heatmap(adata: ad.AnnData, score_cols: list[str], out_dir: Path) -> None:
    if not score_cols:
        return
    data = adata.obs[["tcell_clusters", *score_cols]].copy()
    cluster_means = data.groupby("tcell_clusters")[score_cols].mean().sort_index()
    plt.figure(figsize=(8, max(3, cluster_means.shape[0] * 0.35)))
    sns.heatmap(cluster_means, cmap="coolwarm", center=0, cbar_kws={"label": "Average score"})
    plt.title("T-cell signature scores by cluster")
    plt.xlabel("Signature")
    plt.ylabel("Cluster")
    plt.tight_layout()
    plt.savefig(out_dir / "tcell_signature_heatmap.png", dpi=300)
    plt.close()

def plot_cluster_composition(adata: ad.AnnData, groupby: str, out_dir: Path) -> None:
    if groupby not in adata.obs:
        return
    counts = adata.obs.groupby([groupby, "tcell_clusters"]).size().unstack(fill_value=0)
    if counts.empty:
        return
    fractions = counts.div(counts.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    if fractions.empty:
        return
    fig, ax = plt.subplots(figsize=(max(6, fractions.shape[0] * 0.35), 4))
    fractions.sort_index().plot(kind="bar", stacked=True, ax=ax, colormap="tab20")
    ax.set_ylabel("Fraction of cells")
    ax.set_xlabel(groupby.capitalize())
    ax.set_title(f"T-cell cluster composition by {groupby}")
    ax.legend(title="Cluster", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    plt.tight_layout()
    plt.savefig(out_dir / f"tcell_cluster_proportions_by_{groupby}.png", dpi=300)
    plt.close(fig)

def rank_cluster_markers(adata: ad.AnnData, out_dir: Path, top_n: int = 5, n_genes: int = 50) -> None:
    if "tcell_clusters" not in adata.obs:
        return
    sc.tl.rank_genes_groups(adata, groupby="tcell_clusters", method="wilcoxon", n_genes=n_genes, use_raw=False)
    markers = sc.get.rank_genes_groups_df(adata, None)
    markers.to_csv(out_dir / "tcell_cluster_markers.csv", index=False)

    cluster_series = adata.obs["tcell_clusters"]
    if pd.api.types.is_categorical_dtype(cluster_series):
        categories = list(cluster_series.cat.categories)
    else:
        categories = sorted(cluster_series.unique())

    var_map: OrderedDict[str, list[str]] = OrderedDict()
    for cat in categories:
        genes = markers[markers["group"] == cat]["names"].head(top_n).tolist()
        if genes:
            var_map[str(cat)] = genes

    if var_map:
        sc.pl.dotplot(adata, var_names=var_map, groupby="tcell_clusters", standard_scale="var", save="_tcell_marker_dotplot.png", show=False)

def run_step11_tcells() -> None:
    out_dir = config.DEFAULT_OUTPUT_SUBDIRS["t_cells"]
    ensure_dirs([out_dir])
    sc.settings.figdir = str(out_dir)

    adata = load_celltype_data()
    t_cells = subset_t_cells(adata)
    preprocess(t_cells)
    score_cols = score_signatures(t_cells)

    cluster_signature_heatmap(t_cells, score_cols, out_dir)
    for group in ["sample", "patient", "tissue", "orig.ident"]:
        plot_cluster_composition(t_cells, group, out_dir)
    rank_cluster_markers(t_cells, out_dir)

    sc.pl.tsne(t_cells, color=["tcell_label", "tcell_clusters"], save="_tcell_labels.png", show=False)
    sc.pl.umap(t_cells, color=score_cols, save="_signature_scores.png", show=False)

    t_cells.write_h5ad(out_dir / "t_cells_annotated.h5ad")
    LOGGER.info("Saved T-cell annotations to %s", out_dir / "t_cells_annotated.h5ad")
