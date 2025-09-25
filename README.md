# Gastric-Cancer-scRNA-seq-Analysis-Pipeline---GSE163558

# Gastric Cancer scRNA-seq Analysis Pipeline

This project is a Python reimplementation of the gastric cancer single-cell RNA sequencing analysis pipeline based on the GSE163558 dataset, including TCGA data integration.

## Project Overview

This is a comprehensive single-cell RNA sequencing analysis pipeline for studying the tumor microenvironment in gastric cancer (STAD - Stomach Adenocarcinoma). The project includes a complete workflow from raw data processing to advanced analysis, encompassing:

- Data quality control and preprocessing
- Batch effect correction (Harmony)
- Cell type annotation
- Malignant cell scoring
- Gene signature analysis
- Pseudotime analysis
- Survival analysis
- CNV (Copy Number Variation) analysis
- T cell subtype annotation
- Cell-cell communication analysis

## Datasets

- **scRNA-seq data**: GSE163558 (containing 10 samples: PT1-3, NT1, LN1-2, O1, P1, Li1-2)
- **TCGA data**: TCGA-STAD (gastric cancer bulk RNA-seq data for validation and survival analysis)

## Project Structure

```
├── GSE163558_RAW/           # Raw 10X Genomics data
├── input/                    # TCGA input data
├── python_pipeline/          # Python analysis pipeline
│   ├── scripts/             # Analysis step scripts
│   ├── config.py            # Configuration parameters
│   ├── preprocessing.py     # Data preprocessing
│   ├── clustering.py        # Clustering analysis
│   ├── malignant_scoring.py # Malignant cell scoring
│   └── ...
├── python_outputs/           # Analysis results output
├── scRNA_scripts/            # R language auxiliary scripts
├── cache/                    # Cache files and temporary data
└── step*.R                   # R language analysis steps
```

## Environment Setup

### System Requirements

- Python 3.8+
- R 4.0+
- Conda environment management recommended

### Python Dependencies

Main dependencies are listed in `python_pipeline/requirements.txt`:

```bash
pip install -r python_pipeline/requirements.txt
```

### R Dependencies

Run `install_packages.R` to install required R packages.

## Usage

### 1. Data Preparation

Ensure raw data is located in the `GSE163558_RAW/` directory and TCGA data in the `input/` directory.

### 2. Run Complete Analysis Pipeline

The project contains 12 main analysis steps that can be run sequentially:

```bash
# Navigate to Python pipeline directory
cd python_pipeline

# Step 1-2: Quality control and batch effect correction
python scripts/step01_qc_harmony.py

# Step 3: Cell type annotation
python scripts/step03_celltype.py

# Step 4: Feature visualization
python scripts/step04_featureplots.py

# Step 5: Cell proportion statistics
python scripts/step05_proportions.py

# Step 6: TCGA integration and malignant cell identification
python scripts/step06_tcga_malignant.py

# Step 7: Gene signature analysis
python scripts/step07_signatures.py

# Step 8: Pseudotime analysis
python scripts/step08_pseudotime.py

# Step 9: Survival analysis
python scripts/step09_survival.py

# Step 10: CNV analysis
python scripts/step10_cnv.py

# Step 11: T cell subtype analysis
python scripts/step11_tcells.py

# Step 12: Cell communication analysis
python scripts/step12_cell_communication.py
```

### 3. Output Results

All analysis results will be saved in the `python_outputs/` directory, organized by step numbers:

- `1-QC/`: Quality control results
- `2-harmony/`: Batch-corrected data
- `3-celltype/`: Cell type annotation results
- `4-plot/`: Visualization plots
- `5-prop/`: Cell proportion statistics
- `6-malignant/`: Malignant cell analysis
- `7-signatures/`: Gene signature analysis
- `8-pseudotime/`: Pseudotime analysis
- `9-survival/`: Survival analysis
- `10-cnv/`: CNV analysis
- `11-tcells/`: T cell subtype analysis
- `12-cellchat/`: Cell communication analysis

## Analysis Workflow Details

### 1. Quality Control and Preprocessing
- Load 10X Genomics data
- Filter low-quality cells (genes < 300)
- Filter mitochondrial/ribosomal gene proportions
- Data normalization and HVG selection
- Harmony batch effect correction

### 2. Cell Type Annotation
- Initial annotation based on known marker genes
- Precise cell type identification using SingleR algorithm

### 3. Malignant Cell Identification
- Comparison with TCGA bulk data
- Malignant cell identification based on gene expression profiles
- Calculation of malignancy scores

### 4. Advanced Analysis
- **Signatures**: Immune, stem cell, cell cycle functional analysis
- **Pseudotime**: Cell differentiation trajectory reconstruction
- **Survival**: Prognosis analysis based on cell proportions and signatures
- **CNV**: Copy number variation detection
- **T cell subtypes**: Detailed CD4+/CD8+ T cell classification
- **Cell communication**: Ligand-receptor interaction analysis

## Key Output Files

- `*.h5ad`: Single-cell data in AnnData format (scanpy standard)
- `*.csv`: Statistical results and marker gene lists
- `*.png/*.pdf`: Visualization plots and heatmaps
- `timings/*.jsonl`: Runtime records for each step

## References and Acknowledgments

This project is based on the following public datasets and methods:

- GSE163558: Gastric cancer single-cell RNA sequencing data
- TCGA-STAD: Gastric cancer bulk RNA sequencing data
- Seurat/Signac: Original R implementation
- Scanpy: Python single-cell analysis framework
- Harmony: Batch effect correction
- InferCNV: CNV analysis

## License

This project is for academic research purposes only. Please comply with relevant data usage agreements.

## Contact

Please submit issues on GitHub if you have any questions
