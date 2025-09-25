"""Gene marker collections translated from the original R scripts."""

from __future__ import annotations

from typing import Dict, List

BASIC_MARKER_PANEL: List[str] = [
    "EPCAM", "KRT19", "CLDN4",
    "PECAM1", "COL1A2", "VWF",
    "CD3D", "CD3E", "CD8A", "CD4", "CD2",
    "LUM", "FGF7", "MME",
    "AIF1", "C1QC", "C1QB", "LYZ",
    "MKI67", "STMN1", "PCNA",
    "CPA3", "CST3", "KIT", "TPSAB1", "TPSB2",
    "GOS2", "S100A9", "S100A8", "CXCL8",
    "KLRD1", "GNLY", "KLRF1", "AREG", "XCL2", "HSPA6",
    "MS4A1", "CD19", "CD79A", "IGHG1", "MZB1", "SDC1",
    "CSF1R", "CSF3R", "CD68",
]

CELLTYPE_MAPPING = {
    0: "T",
    1: "T",
    2: "NK",
    3: "B",
    4: "Myeloid",
    5: "Epithelial",
    6: "Myeloid",
    7: "B",
    8: "Fibro",
    9: "Proliferative",
    10: "Endothelial",
    11: "Epithelial",
    12: "Epithelial",
    13: "Myeloid",
}

MALIGNANT_SIGNATURES = {
    "tumor_chars": [
        "CLDN7", "CLDN4", "TFF3", "LIPF", "MUC5AC", "GKN1", "PGC",
    ],
    "stemness": ["CD44", "PROM1", "LGR5", "SOX2", "TFRC", "CXCR4", "JAG1"],
}

PROLIFERATION_GENESET = ["MKI67", "IGF1", "ITGB2", "PDGFC", "JAG1", "PHGDH"]
MIGRATION_GENESET = ["VIM", "SNAI1", "MMP9", "AREG", "ARID5B", "FAT1"]

T_CELL_SIGNATURES: Dict[str, List[str]] = {
    "CD4": ["CD3D", "CD3E", "CD4"],
    "CD8": ["CD3D", "CD3E", "CD8A", "CD8B"],
    "Treg": ["FOXP3", "IL2RA", "CTLA4", "IKZF2"],
    "Exhausted": ["LAG3", "PDCD1", "TIGIT", "HAVCR2"],
    "Cytotoxic": ["PRF1", "GZMB", "NKG7", "GNLY"],
    "Memory": ["CCR7", "SELL", "TCF7", "LEF1"],
    "MAIT": ["SLC4A10", "ZBTB16", "RORA"],
}

