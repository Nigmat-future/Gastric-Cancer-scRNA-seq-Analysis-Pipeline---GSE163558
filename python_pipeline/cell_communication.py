"""Cell-cell communication analysis using Squidpy."""

from __future__ import annotations

from pathlib import Path
import os

import anndata as ad
import pandas as pd
try:
    import squidpy as sq  # type: ignore
except Exception:  # pragma: no cover - optional dependency handling
    sq = None  # type: ignore

from . import config
from .utils.io import ensure_dirs
from .utils.logging import get_logger

LOGGER = get_logger(__name__)


def load_tcell_data() -> ad.AnnData:
    path = config.DEFAULT_OUTPUT_SUBDIRS["t_cells"] / "t_cells_annotated.h5ad"
    if not path.exists():
        raise FileNotFoundError("Run step11_tcells.py first to create t_cells_annotated.h5ad")
    return ad.read_h5ad(path)


def load_full_celltype_data() -> ad.AnnData:
    """Load full annotated dataset across all cell types.

    Tries common export locations from the celltype step, falls back to T cell data.
    """
    candidates: list[Path] = []
    try:
        # Likely output path from step3 celltype
        candidates.append(config.DEFAULT_OUTPUT_SUBDIRS["celltype"] / "celltype_annotated.h5ad")
        candidates.append(config.DEFAULT_OUTPUT_SUBDIRS["celltype"] / "annotated.h5ad")
        candidates.append(config.DEFAULT_OUTPUT_SUBDIRS["celltype"] / "celltype.h5ad")
    except Exception:
        pass

    # Try signatures dir as a very last resort
    try:
        candidates.append(config.DEFAULT_OUTPUT_SUBDIRS["signatures"] / "celltype_annotated.h5ad")
    except Exception:
        pass

    for p in candidates:
        if p.exists():
            LOGGER.info("Loading all-cell annotated data from %s", p)
            return ad.read_h5ad(p)

    LOGGER.warning("Full celltype annotated data not found; falling back to T-cell subset")
    return load_tcell_data()


def _pick_group_key(adata: ad.AnnData) -> str:
    """Select a categorical obs key for grouping ligand-receptor analysis."""
    preferred = ["celltype", "malignant_clusters", "tcell_label"]
    for key in preferred:
        if key in adata.obs and pd.api.types.is_categorical_dtype(adata.obs[key]):
            return key
    # fallback: first categorical column
    for col in adata.obs.columns:
        if pd.api.types.is_categorical_dtype(adata.obs[col]):
            return col
    # as a last resort, coerce leiden to category if present
    for key in ["leiden", "clusters", "louvain"]:
        if key in adata.obs:
            adata.obs[key] = adata.obs[key].astype("category")
            return key
    raise RuntimeError("No categorical group key found in adata.obs for ligand-receptor analysis")


def run_cell_communication() -> None:
    out_dir = config.DEFAULT_OUTPUT_SUBDIRS["cellchat"]
    ensure_dirs([out_dir])

    # Prefer full dataset across all cell types
    adata = load_full_celltype_data()
    group_key = _pick_group_key(adata)

    if sq is None:
        msg = (
            "Squidpy is not installed. Install with either 'conda install -c conda-forge squidpy' "
            "or 'pip install squidpy'. If internet is unavailable, this step cannot run."
        )
        LOGGER.error(msg)
        raise ModuleNotFoundError(msg)

    # Determine compute settings
    n_perms = 5000
    n_jobs_env = os.getenv("LR_N_JOBS")
    try:
        n_jobs = int(n_jobs_env) if n_jobs_env is not None else min(os.cpu_count() or 1, 8)
    except Exception:
        n_jobs = min(os.cpu_count() or 1, 8)

    LOGGER.info("Running Squidpy ligand-receptor analysis (OmniPath) with group_key=%s, n_perms=%d, n_jobs=%d", group_key, n_perms, n_jobs)

    # Optional progress bar via tqdm_joblib if joblib backend is used under the hood
    used_progress = False
    try:
        from tqdm.auto import tqdm  # type: ignore
        from tqdm_joblib import tqdm_joblib  # type: ignore
        with tqdm_joblib(tqdm(total=n_perms, desc="LigRec permutations", unit="perm")):
            try:
                sq.gr.ligrec(adata, n_perms=n_perms, use_raw=False, cluster_key=group_key, n_jobs=n_jobs)
            except TypeError:
                sq.gr.ligrec(adata, n_perms=n_perms, use_raw=False, cluster_key=group_key)
        used_progress = True
    except Exception:
        pass

    if not used_progress:
        # Fallback execution without progress bar
        try:
            sq.gr.ligrec(adata, n_perms=n_perms, use_raw=False, cluster_key=group_key, n_jobs=n_jobs)
        except TypeError:
            sq.gr.ligrec(adata, n_perms=n_perms, use_raw=False, cluster_key=group_key)

    # Save annotated results
    adata.write_h5ad(out_dir / "ligrec_allcells.h5ad")

    # Try to export summary tables if available
    try:
        lr_uns = adata.uns.get("ligrec", {})
        # Common keys seen in Squidpy versions
        for key_name in ["pvalues", "means", "results", "summary"]:
            if key_name in lr_uns:
                df = lr_uns[key_name]
                if hasattr(df, "to_csv"):
                    df.to_csv(out_dir / f"ligrec_{key_name}.csv")
    except Exception:
        LOGGER.warning("Could not export ligand-receptor tables; continuing with plots")

    # Heatmap with FDR/p-value threshold where supported via alpha
    try:
        sq.pl.ligrec(
            adata,
            source_groups=group_key,
            target_groups=group_key,
            alpha=0.05,
            save=out_dir / "ligrec_heatmap_allcells.pdf",
            show=False,
        )
    except Exception:
        LOGGER.warning("Ligrec heatmap plotting failed")

    # If both malignant and T labels present, plot a focused cross-talk heatmap
    try:
        has_malignant = "malignant_clusters" in adata.obs
        has_tcell = "tcell_label" in adata.obs
        if has_malignant and has_tcell:
            sq.pl.ligrec(
                adata,
                source_groups="malignant_clusters",
                target_groups="tcell_label",
                alpha=0.05,
                save=out_dir / "ligrec_heatmap_malignant_to_tcell.pdf",
                show=False,
            )
    except Exception:
        LOGGER.warning("Cross-group ligrec plotting failed")

