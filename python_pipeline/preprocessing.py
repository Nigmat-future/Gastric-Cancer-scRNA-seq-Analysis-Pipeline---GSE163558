"""Preprocessing helpers for single-cell data (QC, integration, Harmony)."""

from __future__ import annotations

import json
import dataclasses
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scanpy.external.pp import harmony_integrate

from . import config
from .utils.io import ensure_dirs
from .utils.logging import get_logger, time_block

LOGGER = get_logger(__name__)


def prepare_raw_matrix_structure(raw_dir: Path | None = None) -> List[Path]:
    """Ensure each 10X sample has its own folder with matrix.mtx.gz etc.

    Mimics the R script that reorganises files downloaded from GEO so that
    scanpy can read them directly. Returns the list of created sample folders.
    """

    raw_dir = raw_dir or config.RAW_10X_DIR
    if not raw_dir.exists():
        raise FileNotFoundError(f"Raw directory not found: {raw_dir}")

    folders: List[Path] = []
    for gsm_folder in sorted(raw_dir.glob("GSM*")):
        if gsm_folder.is_dir():
            # Directory already structured.
            folders.append(gsm_folder)
            continue

    if folders:
        LOGGER.info("Detected %d pre-existing 10X folders", len(folders))
        return folders

    gz_files = sorted(raw_dir.glob("GSM*.gz"))
    if not gz_files:
        raise RuntimeError("No GSM gz files found to reorganise")

    grouped: Dict[str, List[Path]] = {}
    for path in gz_files:
        stem = path.name.split("_")[0:2]
        key = "_".join(stem)
        grouped.setdefault(key, []).append(path)

    for key, paths in grouped.items():
        sample_dir = raw_dir / key
        sample_dir.mkdir(exist_ok=True)
        mapping = {
            "barcodes": "barcodes.tsv.gz",
            "features": "features.tsv.gz",
            "genes": "features.tsv.gz",
            "matrix": "matrix.mtx.gz",
        }
        for src in paths:
            token = src.name.split("_")[-1].split(".")[0].lower()
            target_name = mapping.get(token)
            if target_name is None:
                LOGGER.warning("Unknown file token for %s", src.name)
                continue
            dst = sample_dir / target_name
            if not dst.exists():
                src.replace(dst)
        folders.append(sample_dir)
        LOGGER.info("Prepared 10X folder %s", sample_dir)
    return folders


@dataclass(slots=True)
class QCThresholds:
    min_genes: int = config.MIN_FEATURES
    min_cells: int = config.MIN_CELLS
    max_mito: float = 25.0
    min_ribo: float = 3.0
    max_hb: float = 1.0


def read_10x_samples(sample_dirs: Iterable[Path]) -> Dict[str, ad.AnnData]:
    datasets: Dict[str, ad.AnnData] = {}
    for sample_dir in sample_dirs:
        sample_id = sample_dir.name
        LOGGER.info("Reading sample %s", sample_id)
        adata = sc.read_10x_mtx(sample_dir, cache=True)
        adata.var_names_make_unique()
        adata.obs["orig.ident"] = sample_id
        datasets[sample_id] = adata
    return datasets


def merge_samples(datasets: Dict[str, ad.AnnData]) -> ad.AnnData:
    LOGGER.info("Merging %d samples", len(datasets))
    adata = ad.concat(
        datasets,
        label="orig.ident",
        join="outer",
        index_unique=None,
    )
    return adata


def annotate_metadata(adata: ad.AnnData) -> None:
    obs = adata.obs.copy()
    tokens = obs["orig.ident"].str.split("_", expand=True)
    obs["group"] = tokens[1]

    series = obs["orig.ident"].astype(str)
    sample_map = {
        "PT": "GC",
        "NT": "adjacent_nontumor",
        "LN": "GC_lymph_metastasis",
        "O1": "GC_ovary_metastasis",
        "P1": "GC_peritoneum_metastasis",
        "Li": "GC_liver_metastasis",
    }
    obs["sample"] = series
    for key, value in sample_map.items():
        obs["sample"] = obs["sample"].str.replace(rf"GSM\\d+_{key}\\d*", value, regex=True)

    patient_map = {
        "GSM5004180_PT1": "Patient1",
        "GSM5004188_Li1": "Patient1",
        "GSM5004181_PT2": "Patient2",
        "GSM5004183_NT1": "Patient2",
        "GSM5004186_O1": "Patient3",
        "GSM5004185_LN2": "Patient5",
        "GSM5004187_P1": "Patient6",
    }
    obs["patient"] = series.map(patient_map).fillna("Patient4")

    obs["tissue"] = np.select(
        [
            obs["orig.ident"].isin([
                "GSM5004180_PT1",
                "GSM5004181_PT2",
                "GSM5004182_PT3",
            ]),
            obs["orig.ident"].isin(["GSM5004183_NT1"]),
        ],
        ["PT", "NT"],
        default="M",
    )
    adata.obs = obs


def compute_qc_metrics(adata: ad.AnnData) -> None:
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=[],
        inplace=True,
        percent_top=None,
        log1p=False,
    )
    mito_mask = adata.var_names.str.match(config.MITO_PATTERN)
    ribo_mask = adata.var_names.str.match(config.RIBO_PATTERN)
    hb_mask = adata.var_names.str.match(config.HB_PATTERN)

    adata.obs["percent_mito"] = (adata[:, mito_mask].X.sum(axis=1).A1 / adata.obs["total_counts"]) * 100
    adata.obs["percent_ribo"] = (adata[:, ribo_mask].X.sum(axis=1).A1 / adata.obs["total_counts"]) * 100
    adata.obs["percent_hb"] = (adata[:, hb_mask].X.sum(axis=1).A1 / adata.obs["total_counts"]) * 100


def filter_cells_and_genes(adata: ad.AnnData, *, thresholds: QCThresholds | None = None) -> ad.AnnData:
    thresholds = thresholds or QCThresholds()
    LOGGER.info("Filtering cells using thresholds: %s", json.dumps(dataclasses.asdict(thresholds)))
    sc.pp.filter_genes(adata, min_cells=thresholds.min_cells)
    sc.pp.filter_cells(adata, min_genes=thresholds.min_genes)
    mask = (
        (adata.obs["percent_mito"] < thresholds.max_mito)
        & (adata.obs["percent_ribo"] > thresholds.min_ribo)
        & (adata.obs["percent_hb"] < thresholds.max_hb)
    )
    filtered = adata[mask].copy()
    LOGGER.info("Retained %d of %d cells after QC", filtered.n_obs, adata.n_obs)
    return filtered


def normalise_and_hvg(adata: ad.AnnData) -> None:
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=config.N_HVG, subset=False)


def scale_and_reduce(adata: ad.AnnData) -> None:
    sc.pp.scale(adata)
    sc.tl.pca(adata, n_comps=config.N_PCS)


def run_harmony(adata: ad.AnnData, batch_key: str = "orig.ident") -> None:
    LOGGER.info("Running Harmony integration on key=%s", batch_key)
    harmony_integrate(adata, key=batch_key)
    sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_pcs=config.N_PCS)
    sc.tl.umap(adata)
    sc.tl.tsne(adata, use_rep="X_pca_harmony", random_state=config.SEED_TSNE)
    for res in config.LEIDEN_RESOLUTIONS:
        sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res}")


def preprocess_pipeline(output_dir: Path | None = None) -> ad.AnnData:
    """Perform the full QC + Harmony workflow analogous to step1_2."""

    output_dir = output_dir or config.OUTPUT_DIR
    ensure_dirs([output_dir, *config.DEFAULT_OUTPUT_SUBDIRS.values()])
    timings_file = config.DEFAULT_OUTPUT_SUBDIRS["timings"] / "step01_qc_harmony.jsonl"

    with time_block("prepare_raw_matrix_structure", write_jsonl=timings_file):
        sample_dirs = prepare_raw_matrix_structure()
    with time_block("read_10x_samples", write_jsonl=timings_file):
        datasets = read_10x_samples(sample_dirs)
    with time_block("merge_samples", write_jsonl=timings_file):
        adata = merge_samples(datasets)
    with time_block("annotate_metadata", write_jsonl=timings_file):
        annotate_metadata(adata)
    with time_block("compute_qc_metrics", write_jsonl=timings_file):
        compute_qc_metrics(adata)
    with time_block("filter_cells_and_genes", write_jsonl=timings_file):
        filtered = filter_cells_and_genes(adata)
    with time_block("normalise_and_hvg", write_jsonl=timings_file):
        normalise_and_hvg(filtered)
    with time_block("scale_and_reduce", write_jsonl=timings_file):
        scale_and_reduce(filtered)
    with time_block("run_harmony", write_jsonl=timings_file):
        run_harmony(filtered)

    filtered.write_h5ad(config.DEFAULT_OUTPUT_SUBDIRS["harmony"] / "sce_all_int.h5ad")
    return filtered

