"""TCGA survival analysis using lifelines."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, Tuple

import pandas as pd
import mygene
import matplotlib.pyplot as plt

from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

from . import config
from .tcga_bulk import _input_paths, _sample_groups, log_cpm
from .utils.io import ensure_dirs
from .utils.logging import get_logger

LOGGER = get_logger(__name__)

GENES_OF_INTEREST = ["DST", "PIK3R1", "CXCR4", "GADD45B"]


def load_survival_inputs() -> Tuple[pd.DataFrame, pd.DataFrame]:
    paths = _input_paths()
    survival = pd.read_table(paths['survival'], index_col='sample')
    clinical = pd.read_table(paths['phenotype'], sep='\t')
    for col in ('submitter_id.samples', 'sample', 'sample_id.samples', 'case_submitter_id.samples'):
        if col in clinical.columns:
            clinical = clinical.set_index(col)
            break
    else:
        raise KeyError('Could not locate a sample identifier column in the phenotype file')
    clinical.index = clinical.index.astype(str)
    return survival, clinical


def _map_expression_to_symbols(log_expr: pd.DataFrame) -> pd.DataFrame:
    base_ids = log_expr.index.to_series().str.replace(r"\.\d+$", '', regex=True)
    mg = mygene.MyGeneInfo()
    try:
        hits = mg.querymany(base_ids.unique().tolist(), scopes='ensembl.gene', fields='symbol', species='human')
    except Exception as exc:
        LOGGER.warning('Gene symbol lookup failed: %s', exc)
        return log_expr

    mapping: Dict[str, str] = {}
    for hit in hits:
        if not isinstance(hit, dict):
            continue
        query = hit.get('query')
        symbol = hit.get('symbol')
        if query and symbol:
            mapping.setdefault(query, symbol.upper())

    mapped = base_ids.map(mapping)
    if mapped.isna().all():
        LOGGER.warning('Unable to map Ensembl IDs to gene symbols; keeping original indices')
        return log_expr

    keep = mapped.notna()
    filtered = log_expr.loc[keep].copy()
    filtered.index = mapped.loc[keep]
    aggregated = filtered.groupby(filtered.index, sort=False).mean()
    return aggregated

def prepare_expression() -> pd.DataFrame:
    paths = _input_paths()
    counts = pd.read_table(paths["counts"], index_col=0)
    groups = _sample_groups(counts.columns)
    tumor_cols = groups[groups == "tumor"].index
    expression = counts[tumor_cols]
    log_expr = log_cpm(expression)
    log_expr = log_expr.loc[:, ~log_expr.columns.duplicated()]
    log_expr.columns = log_expr.columns.str.slice(0, 16)
    log_expr = _map_expression_to_symbols(log_expr)
    return log_expr


def prepare_metadata() -> pd.DataFrame:
    survival, clinical = load_survival_inputs()
    meta = survival.join(clinical, how='left')
    meta = meta.loc[~meta.index.duplicated(keep='first')]

    rename_map = {}
    if 'OS' in meta.columns:
        rename_map['OS'] = 'event'
    if 'OS.time' in meta.columns:
        rename_map['OS.time'] = 'time'
    if 'pathologic_M' in meta.columns:
        rename_map['pathologic_M'] = 'M'
    elif 'ajcc_pathologic_m.diagnoses' in meta.columns:
        rename_map['ajcc_pathologic_m.diagnoses'] = 'M'
    if 'pathologic_N' in meta.columns:
        rename_map['pathologic_N'] = 'N'
    elif 'ajcc_pathologic_n.diagnoses' in meta.columns:
        rename_map['ajcc_pathologic_n.diagnoses'] = 'N'
    if 'pathologic_T' in meta.columns:
        rename_map['pathologic_T'] = 'T'
    elif 'ajcc_pathologic_t.diagnoses' in meta.columns:
        rename_map['ajcc_pathologic_t.diagnoses'] = 'T'

    meta = meta.rename(columns=rename_map)

    if 'event' not in meta.columns or 'time' not in meta.columns:
        missing = {'event', 'time'} - set(meta.columns)
        raise KeyError(f'Missing survival columns: {missing}')

    keep_cols = ['event', 'time']
    for col in ('M', 'N', 'T'):
        if col in meta.columns:
            keep_cols.append(col)
    meta = meta[keep_cols]

    meta['event'] = pd.to_numeric(meta['event'], errors='coerce')
    meta['time'] = pd.to_numeric(meta['time'], errors='coerce')
    meta = meta.dropna(subset=['event', 'time'])
    meta['time'] = meta['time'] / 30.0
    meta.index = meta.index.str.slice(0, 16)
    return meta


def align_expression_metadata(expr: pd.DataFrame, meta: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    shared = expr.columns.intersection(meta.index)
    expr = expr.loc[:, shared]
    meta = meta.loc[shared]
    return expr, meta


def gene_high_low(expr: pd.Series) -> pd.Series:
    q75 = expr.quantile(0.75)
    return pd.Series(pd.cut(expr, bins=[-float("inf"), q75, float("inf")], labels=["low", "high"]), index=expr.index)


def run_survival_analysis(genes: Iterable[str] = GENES_OF_INTEREST) -> None:
    out_dir = config.DEFAULT_OUTPUT_SUBDIRS["survival"]
    ensure_dirs([out_dir])

    expr = prepare_expression()
    meta = prepare_metadata()
    expr, meta = align_expression_metadata(expr, meta)

    expr.to_csv(out_dir / "expression_logcpm.csv")
    meta.to_csv(out_dir / "metadata.csv")

    results = []
    for gene in genes:
        if gene not in expr.index:
            LOGGER.warning("Gene %s not found in expression matrix", gene)
            continue
        groups = gene_high_low(expr.loc[gene])
        data = meta.join(groups.rename('group'), how='inner')
        data = data.dropna(subset=['group'])
        if data['group'].nunique() < 2:
            LOGGER.warning('Gene %s does not split into two groups', gene)
            continue
        km = KaplanMeierFitter()
        fig, ax = plt.subplots(figsize=(6, 4))
        for label, color in [('high', '#d62728'), ('low', '#1f77b4')]:
            subset = data[data['group'] == label]
            if subset.empty:
                continue
            km.fit(subset['time'], event_observed=subset['event'], label=label.capitalize())
            km.plot(ax=ax, ci_show=False, color=color)
        ax.set_title(f'{gene} high vs low')
        ax.set_xlabel('Time (months)')
        ax.set_ylabel('Survival probability')
        ax.legend(frameon=False)
        ax.grid(alpha=0.3, linestyle='--', linewidth=0.5)
        fig.tight_layout()
        fig_path_pdf = out_dir / f'km_{gene}.pdf'
        fig_path_png = out_dir / f'km_{gene}.png'
        fig.savefig(fig_path_pdf)
        fig.savefig(fig_path_png, dpi=300)
        plt.close(fig)

        lr = logrank_test(
            data[data['group'] == 'high']['time'],
            data[data['group'] == 'low']['time'],
            event_observed_A=data[data['group'] == 'high']['event'],
            event_observed_B=data[data['group'] == 'low']['event'],
        )
        results.append({'gene': gene, 'p_value': lr.p_value, 'test_statistic': lr.test_statistic})

    if results:
        pd.DataFrame(results).to_csv(out_dir / "survival_tests.csv", index=False)

