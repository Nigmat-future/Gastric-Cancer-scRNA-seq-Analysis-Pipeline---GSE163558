"""CNV inference using infercnvpy as a Python analogue of CopyKAT."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import infercnvpy
import numpy as np
import pandas as pd
from pybiomart import Dataset

from . import config
from .utils.io import ensure_dirs
from .utils.logging import get_logger

LOGGER = get_logger(__name__)


def load_epithelial_data() -> ad.AnnData:
    path = config.DEFAULT_OUTPUT_SUBDIRS["malignant"] / "tumor_epithelial.h5ad"
    if not path.exists():
        raise FileNotFoundError("Run step06_tcga_malignant.py first to create tumor_epithelial.h5ad")
    return ad.read_h5ad(path)


def annotate_gene_positions(adata: ad.AnnData) -> ad.AnnData:
    required = {'chromosome', 'start', 'end'}
    if required.issubset(adata.var.columns):
        return adata

    dataset = Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')
    var = adata.var.copy()
    var['symbol'] = var.index.astype(str)
    var['symbol_upper'] = var['symbol'].str.upper()
    var['ensembl_id'] = var['symbol'].str.split('.').str[0]

    def _query(filter_name: str, values: list[str]) -> pd.DataFrame:
        if not values:
            return pd.DataFrame()
        results = []
        step = 200
        for i in range(0, len(values), step):
            subset = values[i:i + step]
            try:
                df = dataset.query(
                    attributes=['ensembl_gene_id', 'hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'],
                    filters={filter_name: subset},
                )
            except Exception as exc:
                LOGGER.warning('Failed to fetch gene annotations from Ensembl (%s batch): %s', filter_name, exc)
                return pd.DataFrame()
            results.append(df)
        return pd.concat(results, ignore_index=True) if results else pd.DataFrame()

    symbol_df = _query('hgnc_symbol', var['symbol_upper'].dropna().unique().tolist())
    ensembl_candidates = var.loc[var['ensembl_id'].str.startswith('ENSG', na=False), 'ensembl_id'].unique().tolist()
    ensembl_df = _query('ensembl_gene_id', ensembl_candidates)
    lookup = pd.concat([symbol_df, ensembl_df], ignore_index=True) if not symbol_df.empty or not ensembl_df.empty else pd.DataFrame()
    if lookup.empty:
        LOGGER.warning('Unable to retrieve gene positions from Ensembl; skipping CNV inference')
        return adata[:, []].copy()

    lookup['chromosome_name'] = lookup['chromosome_name'].astype(str)
    valid_chr = {str(i) for i in range(1, 23)} | {'X', 'Y'}
    lookup = lookup[lookup['chromosome_name'].isin(valid_chr)]

    chrom = []
    start = []
    end = []
    for gene, upper, ens in zip(var['symbol'], var['symbol_upper'], var['ensembl_id']):
        row = lookup[lookup['hgnc_symbol'].str.upper() == upper]
        if row.empty and isinstance(ens, str):
            row = lookup[lookup['ensembl_gene_id'] == ens]
        if row.empty:
            chrom.append(np.nan)
            start.append(np.nan)
            end.append(np.nan)
        else:
            rec = row.iloc[0]
            chrom.append(rec['chromosome_name'])
            start.append(int(rec['start_position']))
            end.append(int(rec['end_position']))

    adata.var['chromosome'] = chrom
    adata.var['start'] = start
    adata.var['end'] = end

    mask = adata.var['chromosome'].notna()
    dropped = (~mask).sum()
    if dropped:
        LOGGER.warning('Dropping %d genes without genomic coordinates for CNV inference', int(dropped))
    return adata[:, mask].copy()



def run_cnv_inference() -> None:
    out_dir = config.DEFAULT_OUTPUT_SUBDIRS["cnv"]
    ensure_dirs([out_dir])

    adata = load_epithelial_data()
    if "malignant_label" not in adata.obs:
        raise KeyError("Expected malignant_label column in AnnData.obs")

    adata = annotate_gene_positions(adata)
    if adata.n_vars == 0:
        LOGGER.warning('Skipping CNV inference because no genomic annotations were available')
        return

    LOGGER.info("Running infercnvpy using nonmalignant epithelial cells as reference")
    infercnvpy.tl.infercnv(
        adata,
        reference_key="malignant_label",
        reference_cat=["nonmalignant"],
        inplace=True,
    )

    adata.write_h5ad(out_dir / "infercnv.h5ad")
    infercnvpy.pl.chromosome_heatmap(
        adata,
        groupby="malignant_label",
        save=out_dir / "cnv_heatmap.pdf",
        show=False,
    )

