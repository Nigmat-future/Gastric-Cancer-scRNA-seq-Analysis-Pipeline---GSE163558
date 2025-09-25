"""Statistical helpers mirroring bespoke R utilities."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List

import numpy as np
import pandas as pd
from scipy import stats


@dataclass(slots=True)
class GeneTestResult:
    gene: str
    pvalue: float
    statistic: float
    comparison: str

    @property
    def significance_label(self) -> str:
        if self.pvalue > 0.05:
            return "ns"
        if self.pvalue > 0.01:
            return "*"
        if self.pvalue > 0.001:
            return "**"
        return "****"


def two_group_ttests(
    adata,
    genes: Iterable[str],
    groupby: str,
    groups: Iterable[str],
    *,
    require_positive: bool = False,
) -> List[GeneTestResult]:
    """Replicates the singlecell_gene_test function from R."""

    g1, g2 = groups
    mask1 = adata.obs[groupby] == g1
    mask2 = adata.obs[groupby] == g2
    matrix = adata.to_df()
    results: List[GeneTestResult] = []
    for gene in genes:
        if gene not in matrix.columns:
            continue
        expr1 = matrix.loc[mask1, gene].to_numpy()
        expr2 = matrix.loc[mask2, gene].to_numpy()
        if require_positive:
            expr1 = expr1[expr1 > 0]
            expr2 = expr2[expr2 > 0]
        if len(expr1) == 0 or len(expr2) == 0:
            continue
        stat, pval = stats.ttest_ind(expr1, expr2, equal_var=False)
        results.append(GeneTestResult(gene=gene, pvalue=pval, statistic=stat, comparison=f"{g1}_vs_{g2}"))
    return results

