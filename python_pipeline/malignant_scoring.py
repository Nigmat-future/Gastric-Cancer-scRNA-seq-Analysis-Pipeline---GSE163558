"""Malignant epithelial scoring based on TCGA-derived signatures."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import mygene

from . import config
from .statistics import two_group_ttests
from .utils.io import ensure_dirs
from .utils.logging import get_logger
from .utils.markers import MALIGNANT_SIGNATURES, PROLIFERATION_GENESET, MIGRATION_GENESET

LOGGER = get_logger(__name__)

SCORE_THRESHOLD = -0.02


def _gene_set_path() -> Path:
    return config.DEFAULT_OUTPUT_SUBDIRS["tcga"] / "malignant_gene_sets.json"


def load_gene_sets() -> Dict[str, Iterable[str]]:
    path = _gene_set_path()
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    LOGGER.warning("Gene set file %s not found; using built-in fallback", path)
    return {"up": MALIGNANT_SIGNATURES["tumor_chars"], "down": MALIGNANT_SIGNATURES["stemness"]}


def harmonize_gene_sets(gene_sets: Dict[str, Iterable[str]], var_names: Iterable[str]) -> Dict[str, List[str]]:
    var_lookup = {name.upper(): name for name in var_names}
    mg = mygene.MyGeneInfo()
    resolved: Dict[str, List[str]] = {}

    for key, genes in gene_sets.items():
        gene_list = list(genes)
        matched: List[str] = []
        pending: List[str] = []

        for gene in gene_list:
            candidates = [gene, gene.split('.')[0]]
            found = None
            for candidate in candidates:
                name = var_lookup.get(candidate.upper())
                if name:
                    found = name
                    break
            if found:
                matched.append(found)
            else:
                pending.append(gene)

        pending_map = {gene: gene.split('.')[0] for gene in pending}
        symbol_map: Dict[str, str] = {}
        if pending_map:
            try:
                queries = list(set(pending_map.values()))
                hits = mg.querymany(queries, scopes='ensembl.gene', fields='symbol', species='human')
            except Exception as exc:
                LOGGER.warning('Gene symbol lookup failed for %s: %s', key, exc)
                hits = []
            for hit in hits:
                query = hit.get('query')
                symbol = hit.get('symbol')
                if not query or not symbol:
                    continue
                symbol_map.setdefault(query, symbol)

        for orig, stripped in pending_map.items():
            symbol = symbol_map.get(stripped)
            if symbol:
                name = var_lookup.get(symbol.upper())
                if name:
                    matched.append(name)
                else:
                    LOGGER.debug('Mapped symbol %s for %s not in var_names', symbol, orig)
            else:
                LOGGER.debug('No mapping found for %s', orig)

        deduped = list(dict.fromkeys(matched))
        if len(deduped) < len(gene_list):
            LOGGER.warning('Retained %d/%d genes for %s after harmonization', len(deduped), len(gene_list), key)
        else:
            LOGGER.info('Gene set %s fully matched (%d genes)', key, len(deduped))
        resolved[key] = deduped

    return resolved



def load_celltype_data() -> ad.AnnData:
    path = config.DEFAULT_OUTPUT_SUBDIRS["celltype"] / "sce_celltype.h5ad"
    if not path.exists():
        raise FileNotFoundError("Run step03_celltype.py first to create sce_celltype.h5ad")
    return ad.read_h5ad(path)


def subset_epithelial(adata: ad.AnnData) -> ad.AnnData:
    subset = adata[adata.obs["celltype"] == "Epithelial"].copy()
    LOGGER.info("Selected %d epithelial cells", subset.n_obs)
    return subset


def preprocess_subset(adata: ad.AnnData) -> None:
    # Prepare epithelial cells for malignant scoring with safe normalization.
    use_raw = False
    if adata.raw is not None:
        LOGGER.info('Using .raw counts for epithelial preprocessing')
        adata.X = adata.raw.to_adata().X.copy()
        use_raw = True
    elif 'counts' in adata.layers:
        LOGGER.info("Using 'counts' layer for epithelial preprocessing")
        adata.X = adata.layers['counts'].copy()
        use_raw = True

    adata.var_names_make_unique()
    try:
        has_negative = float(adata.X.min()) < 0
    except Exception:
        has_negative = False
    already_logged = bool(adata.uns.get('log1p')) or has_negative
    LOGGER.info('Normalization guard: use_raw=%s already_logged=%s has_negative=%s', use_raw, already_logged, has_negative)

    if use_raw or not already_logged:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    else:
        LOGGER.info('Input already log-transformed; skipping normalize_total/log1p')

    hv_mask = None
    if has_negative and 'highly_variable' in adata.var:
        mask = adata.var['highly_variable'].to_numpy().astype(bool)
        if mask.sum() == 0:
            LOGGER.warning('Existing highly_variable mask empty; recomputing with pearson residuals')
            sc.pp.highly_variable_genes(adata, n_top_genes=config.N_HVG, subset=False, flavor='pearson_residuals')
            hv_mask = adata.var['highly_variable'].to_numpy().astype(bool)
        else:
            LOGGER.info('Using existing highly_variable mask with %d genes', int(mask.sum()))
            hv_mask = mask
    else:
        sc.pp.highly_variable_genes(adata, n_top_genes=config.N_HVG, subset=False, flavor='seurat_v3')
        hv_mask = adata.var['highly_variable'].to_numpy().astype(bool)

    if hv_mask is None or hv_mask.sum() == 0:
        LOGGER.warning('No highly variable genes detected; using all genes for embeddings')
        hv_mask = np.ones(adata.n_vars, dtype=bool)

    adata.var['highly_variable'] = hv_mask
    sc.pp.scale(adata)
    if np.isnan(adata.X).any():
        LOGGER.warning('NaNs detected after scaling; replacing with zeros')
        adata.X = np.nan_to_num(adata.X, nan=0.0, posinf=0.0, neginf=0.0)
    sc.tl.pca(adata, n_comps=config.N_PCS, use_highly_variable=True)
    sc.pp.neighbors(adata, n_pcs=config.N_PCS)
    sc.tl.umap(adata)
    sc.tl.tsne(adata, random_state=config.SEED_TSNE)


def score_cells(adata: ad.AnnData, up_genes: Iterable[str], down_genes: Iterable[str]) -> None:
    sc.tl.score_genes(adata, gene_list=list(up_genes), score_name="tumor_score")
    sc.tl.score_genes(adata, gene_list=list(down_genes), score_name="normal_score")
    adata.obs["malignant_score"] = adata.obs["tumor_score"] - adata.obs["normal_score"]
    adata.obs["malignant_label"] = np.where(
        adata.obs["malignant_score"] > SCORE_THRESHOLD, "malignant", "nonmalignant"
    )


def module_scores(adata: ad.AnnData) -> None:
    sc.tl.score_genes(adata, gene_list=MALIGNANT_SIGNATURES["stemness"], score_name="stemness_score")
    sc.tl.score_genes(adata, gene_list=PROLIFERATION_GENESET, score_name="proliferation_score")
    sc.tl.score_genes(adata, gene_list=MIGRATION_GENESET, score_name="migration_score")


def plot_scores(adata: ad.AnnData, out_dir: Path) -> None:
    sns.set(style="whitegrid")
    plt.figure(figsize=(6, 4))
    sns.violinplot(
        data=adata.obs,
        x="malignant_label",
        y="malignant_score",
        palette=["red", "navy"],
    )
    plt.tight_layout()
    plt.savefig(out_dir / "malignant_score_violin.pdf")
    plt.close()


def plot_gene_violins(adata: ad.AnnData, genes: Iterable[str], out_dir: Path) -> None:
    df = adata[:, list(genes)].to_df()
    df["label"] = adata.obs["malignant_label"].values
    melted = df.melt(id_vars="label", var_name="gene", value_name="expression")
    plt.figure(figsize=(14, 6))
    sns.violinplot(
        data=melted,
        x="gene",
        y="expression",
        hue="label",
        split=False,
        palette=["red", "navy"],
    )
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(out_dir / "marker_violin.pdf")
    plt.close()


def run_step6_malignant() -> None:
    out_dir = config.DEFAULT_OUTPUT_SUBDIRS["malignant"]
    ensure_dirs([out_dir])
    sc.settings.figdir = str(out_dir)

    gene_sets = load_gene_sets()
    adata = load_celltype_data()
    gene_sets = harmonize_gene_sets(gene_sets, adata.var_names)
    if not gene_sets.get('up'):
        LOGGER.warning('Falling back to built-in tumor signature for malignant scoring')
        gene_sets['up'] = MALIGNANT_SIGNATURES['tumor_chars']
    if not gene_sets.get('down'):
        LOGGER.warning('Falling back to built-in normal signature for malignant scoring')
        gene_sets['down'] = MALIGNANT_SIGNATURES['stemness']
    up_genes, down_genes = gene_sets.get('up', []), gene_sets.get('down', [])

    epi = subset_epithelial(adata)
    preprocess_subset(epi)
    score_cells(epi, up_genes, down_genes)
    module_scores(epi)

    sc.pl.tsne(epi, color=["malignant_label", "malignant_score"], save="_malignant_overview.png", show=False)
    sc.pl.umap(epi, color=["malignant_label", "stemness_score", "proliferation_score", "migration_score"], save="_module_scores.png", show=False)

    plot_scores(epi, out_dir)
    plot_gene_violins(epi, MALIGNANT_SIGNATURES["tumor_chars"], out_dir)

    tests = two_group_ttests(
        epi,
        genes=MALIGNANT_SIGNATURES["tumor_chars"],
        groupby="malignant_label",
        groups=("malignant", "nonmalignant"),
        require_positive=True,
    )
    df = pd.DataFrame(
        {
            "gene": [t.gene for t in tests],
            "pvalue": [t.pvalue for t in tests],
            "statistic": [t.statistic for t in tests],
            "significance": [t.significance_label for t in tests],
        }
    )
    df.to_csv(out_dir / "violin_stats.csv", index=False)

    epi.write_h5ad(out_dir / "tumor_epithelial.h5ad")
    LOGGER.info("Saved malignant-scored AnnData to %s", out_dir / "tumor_epithelial.h5ad")

