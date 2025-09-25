"""Cell composition analysis."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from . import config
from .utils.io import ensure_dirs
from .utils.logging import get_logger

LOGGER = get_logger(__name__)


COLORS = [
    "#3176B7", "#F78000", "#3FA116", "#CE2820", "#9265C1",
    "#885649", "#DD76C5", "#BBBE00", "#41BED1",
]


def load_celltype_data(path: Path | None = None) -> ad.AnnData:
    path = path or (config.DEFAULT_OUTPUT_SUBDIRS["celltype"] / "sce_celltype.h5ad")
    if not path.exists():
        raise FileNotFoundError(f"Annotated AnnData not found: {path}")
    return ad.read_h5ad(path)


def _composition_table(adata: ad.AnnData, groupby: str) -> pd.DataFrame:
    counts = adata.obs.groupby([groupby, "celltype"]).size().reset_index(name="count")
    totals = counts.groupby(groupby)["count"].transform("sum")
    counts["percent"] = counts["count"] / totals
    return counts


def _plot_bar(df: pd.DataFrame, groupby: str, out_file: Path) -> None:
    plt.figure(figsize=(12, 6))
    df_sorted = df.sort_values(groupby)
    sns.barplot(
        data=df_sorted,
        x="percent",
        y=groupby,
        hue="celltype",
        orient="h",
        palette=COLORS,
    )
    plt.xlabel("Relative proportion")
    plt.ylabel("")
    plt.legend(loc="upper center", bbox_to_anchor=(0.5, 1.15), ncol=4)
    plt.tight_layout()
    plt.savefig(out_file)
    plt.close()


def run_step5_proportions() -> None:
    out_dir = config.DEFAULT_OUTPUT_SUBDIRS["proportions"]
    ensure_dirs([out_dir])

    adata = load_celltype_data()
    for group in ["tissue", "sample"]:
        df = _composition_table(adata, group)
        out_path = out_dir / f"composition_{group}.csv"
        df.to_csv(out_path, index=False)
        LOGGER.info("Saved composition table for %s to %s", group, out_path)
        _plot_bar(df, group, out_dir / f"composition_{group}.pdf")

    combined = pd.concat(
        [_composition_table(adata, "tissue"), _composition_table(adata, "sample")],
        keys=["tissue", "sample"],
        names=["level", None],
    )
    combined.to_csv(out_dir / "composition_all_levels.csv")

