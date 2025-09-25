"""Central configuration for the Python reimplementation of the gastric cancer scRNA-seq workflow."""

from __future__ import annotations

from pathlib import Path

# Repository layout ---------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent
RAW_10X_DIR = REPO_ROOT / "GSE163558_RAW"
OUTPUT_DIR = REPO_ROOT / "python_outputs"
REFERENCE_DIR = REPO_ROOT / "reference_data"
TCGA_DIR = OUTPUT_DIR / "tcga"
TIMINGS_DIR = OUTPUT_DIR / "timings"

# Ensure lazily created folders are discoverable without touching the FS at import time.
DEFAULT_OUTPUT_SUBDIRS = {
    "qc": OUTPUT_DIR / "1-QC",
    "harmony": OUTPUT_DIR / "2-harmony",
    "celltype": OUTPUT_DIR / "3-celltype",
    "plots": OUTPUT_DIR / "4-plot",
    "proportions": OUTPUT_DIR / "5-prop",
    "tcga": TCGA_DIR,
    "malignant": OUTPUT_DIR / "6-malignant",
    "signatures": OUTPUT_DIR / "7-signatures",
    "pseudotime": OUTPUT_DIR / "8-pseudotime",
    "survival": OUTPUT_DIR / "9-survival",
    "cnv": OUTPUT_DIR / "10-cnv",
    "t_cells": OUTPUT_DIR / "11-tcells",
    "cellchat": OUTPUT_DIR / "12-cellchat",
    "timings": TIMINGS_DIR,
}

# Analysis defaults --------------------------------------------------------
MIN_FEATURES = 300
MIN_CELLS = 5
MITO_PATTERN = "^MT-"
RIBO_PATTERN = "^RP[SL]"
HB_PATTERN = "^HB(?!P)"

# Dimensionality and clustering settings (mirrors Seurat configuration)
N_HVG = 2000
N_PCS = 30
TSNE_DIMS = 15
UMAP_DIMS = 15
LEIDEN_RESOLUTIONS = [0.1, 0.2, 0.5, 0.8, 1.0]

# Random seeds -------------------------------------------------------------
SEED_GLOBAL = 12345
SEED_TSNE = 321

