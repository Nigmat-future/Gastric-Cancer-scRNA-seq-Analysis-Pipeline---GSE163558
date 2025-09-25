"""TCGA-STAD bulk RNA-seq processing and differential expression."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA

import pandas as pd
from statsmodels.stats.multitest import multipletests
from urllib.request import urlopen, Request
import shutil

from . import config
from .utils.io import ensure_dirs
from .utils.logging import get_logger

LOGGER = get_logger(__name__)


@dataclass(slots=True)
class DEGResult:
    table: pd.DataFrame
    tumor_genes: pd.Index
    normal_genes: pd.Index


def _input_paths(proj: str = "TCGA-STAD") -> Dict[str, Path]:
    base = config.DEFAULT_OUTPUT_SUBDIRS["tcga"] / "input"
    return {
        "counts": base / f"{proj}.htseq_counts.tsv.gz",
        "phenotype": base / f"{proj}.GDC_phenotype.tsv.gz",
        "survival": base / f"{proj}.survival.tsv",
    }


def _download_with_ua(url: str, dest: Path) -> None:
    req = Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urlopen(req, timeout=120) as r, open(dest, "wb") as f:
        shutil.copyfileobj(r, f)


def download_tcga_inputs(proj: str = "TCGA-STAD", *, overwrite: bool = False) -> Dict[str, Path]:
    """Download TCGA inputs from UCSC Xena (gdc-hub) to match expected filenames.

    Files:
    - {proj}.htseq_counts.tsv.gz
    - {proj}.GDC_phenotype.tsv.gz
    - {proj}.survival.tsv
    """
    ensure_dirs([config.DEFAULT_OUTPUT_SUBDIRS["tcga"] / "input"])
    paths = _input_paths(proj)

    base_urls = [
        "https://gdc-hub.s3.us-east-1.amazonaws.com/download",
        "https://gdc.xenahubs.net/download",
        "https://gdc-hub.s3.amazonaws.com/download",
    ]

    for key, path in paths.items():
        if path.exists() and not overwrite:
            LOGGER.info("TCGA %s already exists, skipping: %s", key, path)
            continue
        filename = {
            "counts": f"{proj}.htseq_counts.tsv.gz",
            "phenotype": f"{proj}.GDC_phenotype.tsv.gz",
            "survival": f"{proj}.survival.tsv",
        }[key]
        last_err: Exception | None = None
        for base in base_urls:
            url = f"{base}/{filename}"
            try:
                LOGGER.info("Downloading %s → %s", url, path)
                path.parent.mkdir(parents=True, exist_ok=True)
                _download_with_ua(url, path)
                last_err = None
                break
            except Exception as exc:
                LOGGER.warning("Failed %s: %s", url, exc)
                last_err = exc
        if last_err is not None:
            raise last_err
    return paths


def load_counts_table(path: Path) -> pd.DataFrame:
    return pd.read_table(path, index_col=0)


def _filter_genes(df: pd.DataFrame) -> pd.DataFrame:
    nonzero = df.loc[df.sum(axis=1) > 0]
    expressed = nonzero.loc[(nonzero > 0).sum(axis=1) > 0.5 * nonzero.shape[1]]
    return expressed


def _sample_groups(columns: pd.Index) -> pd.Series:
    codes = columns.str.slice(13, 15).astype(int)
    group = np.where(codes < 10, "tumor", "normal")
    return pd.Series(group, index=columns, name="group")


def log_cpm(df: pd.DataFrame) -> pd.DataFrame:
    counts = df.values
    library_sizes = counts.sum(axis=0)
    cpm = (counts / library_sizes) * 1e6
    log_cpm = np.log2(cpm + 1)
    return pd.DataFrame(log_cpm, index=df.index, columns=df.columns)


def differential_expression(df: pd.DataFrame, groups: pd.Series) -> pd.DataFrame:
    log_df = log_cpm(df)
    tumors = log_df.loc[:, groups == "tumor"]
    normals = log_df.loc[:, groups == "normal"]

    mean_t = tumors.mean(axis=1)
    mean_n = normals.mean(axis=1)
    logfc = mean_t - mean_n

    from scipy import stats

    t_stats, p_values = stats.ttest_ind(
        tumors.T,
        normals.T,
        equal_var=False,
        nan_policy="omit",
    )

    deg = pd.DataFrame(
        {
            "logFC": logfc,
            "t_stat": t_stats,
            "pvalue": p_values,
            "mean_tumor": mean_t,
            "mean_normal": mean_n,
        },
        index=df.index,
    ).dropna()
    deg["FDR"] = multipletests(deg["pvalue"], method="fdr_bh")[1]
    deg.sort_values("logFC", ascending=False, inplace=True)
    return deg




def _plot_pca(log_df: pd.DataFrame, groups: pd.Series, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    pca = PCA(n_components=2, random_state=config.SEED_GLOBAL)
    coords = pca.fit_transform(log_df.T)
    explained = pca.explained_variance_ratio_ * 100
    plot_df = pd.DataFrame(coords, columns=['PC1', 'PC2'], index=log_df.columns)
    plot_df['group'] = groups.loc[plot_df.index].values

    plt.figure(figsize=(6, 5))
    sns.scatterplot(data=plot_df, x='PC1', y='PC2', hue='group', palette={'tumor': '#d62728', 'normal': '#1f77b4'})
    plt.xlabel(f'PC1 ({explained[0]:.1f}% var)')
    plt.ylabel(f'PC2 ({explained[1]:.1f}% var)')
    plt.title('TCGA-STAD PCA (log CPM)')
    plt.tight_layout()
    plt.savefig(out_dir / 'tcga_pca.png', dpi=300)
    plt.close()


def _plot_heatmap(log_df: pd.DataFrame, genes: list[str], groups: pd.Series, out_dir: Path) -> None:
    if not genes:
        return
    out_dir.mkdir(parents=True, exist_ok=True)
    data = log_df.loc[genes].copy()
    data = (data - data.mean(axis=1).values[:, None]) / data.std(axis=1).replace(0, 1).values[:, None]
    lut = {'tumor': '#d62728', 'normal': '#1f77b4'}
    col_colors = groups.map(lut)
    sns.clustermap(
        data,
        row_cluster=True,
        col_cluster=True,
        cmap='RdBu_r',
        figsize=(10, 10),
        col_colors=col_colors,
        xticklabels=False,
        yticklabels=True,
    )
    plt.savefig(out_dir / 'tcga_heatmap.png', dpi=300)
    plt.close()


def _plot_volcano(deg: pd.DataFrame, out_dir: Path, *, fdr_cutoff: float = 0.01, logfc_cutoff: float = 1.0) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    df = deg.copy()
    df['neg_log10_fdr'] = -np.log10(df['FDR'])
    df['status'] = 'non-significant'
    df.loc[(df['FDR'] < fdr_cutoff) & (df['logFC'] > logfc_cutoff), 'status'] = 'up'
    df.loc[(df['FDR'] < fdr_cutoff) & (df['logFC'] < -logfc_cutoff), 'status'] = 'down'
    palette = {'up': '#d62728', 'down': '#1f77b4', 'non-significant': 'lightgrey'}

    plt.figure(figsize=(6, 5))
    sns.scatterplot(data=df, x='logFC', y='neg_log10_fdr', hue='status', palette=palette, edgecolor=None, s=12, linewidth=0)
    plt.axvline(logfc_cutoff, linestyle='--', color='grey', linewidth=1)
    plt.axvline(-logfc_cutoff, linestyle='--', color='grey', linewidth=1)
    plt.axhline(-np.log10(fdr_cutoff), linestyle='--', color='grey', linewidth=1)
    plt.xlabel('log2 Fold Change')
    plt.ylabel('-log10 FDR')
    plt.title('TCGA Tumor vs Normal (Volcano)')
    plt.legend(title='', loc='upper right', frameon=False)
    plt.tight_layout()
    plt.savefig(out_dir / 'tcga_volcano.png', dpi=300)
    plt.close()

def select_top_genes(deg: pd.DataFrame, *, fdr_cutoff: float = 0.01, top_n: int = 50) -> Tuple[pd.Index, pd.Index]:
    up = deg[(deg["logFC"] > 1) & (deg["FDR"] < fdr_cutoff)].head(top_n).index
    down = deg[(deg["logFC"] < -1) & (deg["FDR"] < fdr_cutoff)].head(top_n).index
    return up, down


def run_tcga_pipeline() -> DEGResult:
    ensure_dirs([config.DEFAULT_OUTPUT_SUBDIRS["tcga"] / "input"])
    paths = _input_paths()
    for key, path in paths.items():
        if not path.exists():
            raise FileNotFoundError(f"Required TCGA file missing: {path}")

    counts = load_counts_table(paths["counts"])
    counts_filtered = _filter_genes(counts)
    groups = _sample_groups(counts_filtered.columns)
    log_df = log_cpm(counts_filtered)
    deg = differential_expression(counts_filtered, groups)

    deg_path = config.DEFAULT_OUTPUT_SUBDIRS["tcga"] / "tcga_deg.csv"
    ensure_dirs([deg_path.parent])
    deg.to_csv(deg_path)
    LOGGER.info("Saved TCGA DEG table to %s", deg_path)

    plot_dir = config.DEFAULT_OUTPUT_SUBDIRS["tcga"] / "plots"
    _plot_pca(log_df, groups, plot_dir)

    up, down = select_top_genes(deg)
    genes_for_heatmap = list(up) + list(down)
    _plot_heatmap(log_df, genes_for_heatmap, groups, plot_dir)
    _plot_volcano(deg, plot_dir)

    gene_sets_path = config.DEFAULT_OUTPUT_SUBDIRS["tcga"] / "malignant_gene_sets.json"
    gene_sets_path.write_text(json.dumps({"up": list(up), "down": list(down)}), encoding="utf-8")
    LOGGER.info("Exported malignant gene sets to %s", gene_sets_path)

    return DEGResult(table=deg, tumor_genes=up, normal_genes=down)

