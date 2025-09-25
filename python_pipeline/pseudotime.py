"""Trajectory inference using scanpy DPT to mimic Monocle."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from . import config
from .utils.io import ensure_dirs
from .utils.logging import get_logger

LOGGER = get_logger(__name__)


def load_malignant_clusters() -> ad.AnnData:
    path = config.DEFAULT_OUTPUT_SUBDIRS["signatures"] / "malignant_clusters.h5ad"
    if not path.exists():
        raise FileNotFoundError("Run step07_signatures.py first to create malignant_clusters.h5ad")
    return ad.read_h5ad(path)


def _auto_select_root_index(adata: ad.AnnData) -> int:
    """Choose a root cell index based on an early-state gene signature.

    Strategy:
    - Try to score a stemness/early signature (fallback to epithelial/tumor_chars if missing)
    - Pick the cell with the highest score as the root
    - If scoring fails, fallback to the first cell (index 0)
    """
    # Lazy import to avoid hard dependency if markers module changes
    early_genes: list[str] | None = None
    try:
        from .utils.markers import MALIGNANT_SIGNATURES  # type: ignore

        # Prefer stemness if available, else tumor characteristics
        if isinstance(MALIGNANT_SIGNATURES, dict):
            if "stemness" in MALIGNANT_SIGNATURES:
                early_genes = list(MALIGNANT_SIGNATURES["stemness"])  # type: ignore[index]
            elif "tumor_chars" in MALIGNANT_SIGNATURES:
                early_genes = list(MALIGNANT_SIGNATURES["tumor_chars"])  # type: ignore[index]
    except Exception:
        early_genes = None

    if not early_genes:
        # Minimal safe fallback
        return 0

    # Intersect with available genes
    available = [g for g in early_genes if g in adata.var_names]
    if len(available) == 0:
        return 0

    score_name = "score__auto_early"
    try:
        sc.tl.score_genes(adata, gene_list=available, score_name=score_name)
        scores = adata.obs[score_name].to_numpy()
        root_idx = int(np.nanargmax(scores))
        return root_idx
    except Exception:
        return 0


def _maybe_compute_scvi_latent(adata: ad.AnnData) -> bool:
    """Compute scVI latent representation on GPU if available.

    Returns True if X_scVI was added, else False.
    """
    try:
        import torch  # type: ignore
        import scvi  # type: ignore
    except Exception:
        return False

    if not torch.cuda.is_available():
        return False

    # Determine batch and counts layer if present
    batch_key = None
    for key in ("orig.ident", "sample", "patient", "batch"):
        if key in adata.obs:
            batch_key = key
            break

    layer = "counts" if "counts" in adata.layers else None

    try:
        scvi.settings.seed = getattr(config, "SEED_SCVI", 0)
        scvi.model.SCVI.setup_anndata(adata, layer=layer, batch_key=batch_key)
        n_latent = int(getattr(config, "N_PCS", 30))
        model = scvi.model.SCVI(adata, n_latent=n_latent)
        model.train(max_epochs=int(getattr(config, "SCVI_MAX_EPOCHS", 100)), use_gpu=True)
        adata.obsm["X_scVI"] = model.get_latent_representation()
        LOGGER.info("Computed scVI latent representation on GPU; using X_scVI for neighbors")
        return True
    except Exception:
        LOGGER.warning("scVI training failed; falling back to PCA-based neighbors")
        return False


def run_pseudotime() -> None:
    out_dir = config.DEFAULT_OUTPUT_SUBDIRS['pseudotime']
    ensure_dirs([out_dir])
    sc.settings.figdir = str(out_dir)

    adata = load_malignant_clusters()
    adata.obs_names_make_unique()

    # Neighbors: prefer GPU-accelerated scVI latent if available
    if 'neighbors' not in adata.uns:
        if _maybe_compute_scvi_latent(adata):
            sc.pp.neighbors(adata, use_rep='X_scVI')
        else:
            sc.pp.neighbors(adata, n_pcs=config.N_PCS)
    if 'X_diffmap' not in adata.obsm:
        sc.tl.diffmap(adata)

    # PAGA on malignant clusters if available, otherwise first categorical key
    paga_groups = 'malignant_clusters' if 'malignant_clusters' in adata.obs else None
    if paga_groups is None:
        # find first categorical column
        for col in adata.obs.columns:
            if pd.api.types.is_categorical_dtype(adata.obs[col]):
                paga_groups = col
                break
    if paga_groups is not None:
        sc.tl.paga(adata, groups=paga_groups)
        try:
            sc.pl.paga(adata, color=paga_groups, save='_paga.png', show=False)
        except Exception:
            LOGGER.warning('PAGA plotting failed; continuing')

    # Auto-select root if not provided
    if 'iroot' not in adata.uns:
        try:
            adata.uns['iroot'] = _auto_select_root_index(adata)
        except Exception:
            adata.uns['iroot'] = 0

    # DPT pseudotime
    sc.tl.dpt(adata)

    if 'dpt_pseudotime' not in adata.obs:
        LOGGER.warning('DPT failed to add pseudotime; filling with zeros')
        adata.obs['dpt_pseudotime'] = 0.0

    if 'X_umap' not in adata.obsm:
        try:
            sc.tl.umap(adata, init_pos='paga')
        except Exception:
            sc.tl.umap(adata)
    if 'X_tsne' not in adata.obsm:
        sc.tl.tsne(adata, random_state=config.SEED_TSNE)

    sc.pl.umap(adata, color=['dpt_pseudotime', 'malignant_clusters'], save='_pseudotime.png', show=False)
    sc.pl.tsne(adata, color=['dpt_pseudotime'], save='_pseudotime_tsne.png', show=False)

    # Export PAGA edges (if available)
    try:
        if 'paga' in adata.uns and 'connectivities' in adata.uns['paga']:
            conn = adata.uns['paga']['connectivities']
            if hasattr(conn, 'todense'):
                conn_df = pd.DataFrame(np.asarray(conn.todense()))
            else:
                conn_df = pd.DataFrame(np.asarray(conn))
            # Label rows/cols with group names if available
            if paga_groups is not None and paga_groups in adata.obs:
                cats = list(adata.obs[paga_groups].cat.categories)  # type: ignore[attr-defined]
                conn_df.index = cats
                conn_df.columns = cats
            conn_df.to_csv(out_dir / 'paga_connectivities.csv')
    except Exception:
        LOGGER.warning('Failed to export PAGA connectivities')

    # Export cluster-level DE for branches/groups
    try:
        group_key = paga_groups if paga_groups is not None else 'malignant_clusters'
        if group_key in adata.obs:
            sc.tl.rank_genes_groups(adata, groupby=group_key, method='wilcoxon')
            de_df = sc.get.rank_genes_groups_df(adata, None)
            de_df.to_csv(out_dir / 'branch_rank_genes.csv', index=False)
    except Exception:
        LOGGER.warning('Failed to compute/export branch DE tables')

    adata.write_h5ad(out_dir / 'pseudotime.h5ad')

