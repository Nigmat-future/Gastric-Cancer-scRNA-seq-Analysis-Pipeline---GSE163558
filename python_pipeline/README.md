# Python Single-Cell RNA Sequencing Analysis Pipeline

This is the Python reimplementation of the gastric cancer scRNA-seq analysis pipeline, implementing a complete workflow based on the Scanpy ecosystem.

## Core Dependencies

```txt
# Basic analysis
scanpy, anndata, numpy, pandas, scipy, scikit-learn

# Visualization and statistics
matplotlib, seaborn, session_info
gseapy, statsmodels, lifelines, pingouin

# Advanced analysis
harmonypy (batch correction)
infercnvpy (CNV analysis)
squidpy (cell communication)
mygene (gene ID conversion)
```

## Module Structure

| Module | Description |
|--------|-------------|
| `config.py` | Path configuration and global parameters |
| `preprocessing.py` | 10X data loading, QC, normalization, batch correction |
| `clustering.py` | Dimensionality reduction, clustering, marker gene discovery |
| `annotation.py` | Cell type annotation (SingleR logic) |
| `proportions.py` | Cell composition statistical analysis |
| `tcga_bulk.py` | TCGA-STAD data preprocessing and differential analysis |
| `malignant_scoring.py` | Malignant cell scoring and module scores |
| `marker_analysis.py` | Enrichment analysis, heatmaps, dot plots |
| `signature_analysis.py` | Gene signature analysis |
| `pseudotime.py` | Pseudotime analysis (DPT algorithm) |
| `survival.py` | Survival analysis (Kaplan-Meier, Cox) |
| `cnv.py` | CNV inference analysis |
| `tcell_annotation.py` | T cell subtype annotation |
| `cell_communication.py` | Cell-cell communication analysis |
| `scripts/` | Command-line entry scripts (step01-step12) |
| `utils/` | Shared utility functions |

## Analysis Workflow

### Implemented Steps
1. ✅ **step01_qc_harmony.py** - Quality control and Harmony batch correction
2. ✅ **step03_celltype.py** - Cell type annotation
3. ✅ **step04_featureplots.py** - Feature expression visualization
4. ✅ **step05_proportions.py** - Cell proportion statistics
5. ✅ **step06_tcga_malignant.py** - TCGA integration and malignant cell identification
6. ✅ **step07_signatures.py** - Gene signature analysis
7. ✅ **step08_pseudotime.py** - Pseudotime analysis
8. ✅ **step09_survival.py** - Survival analysis
9. ✅ **step10_cnv.py** - CNV analysis
10. ✅ **step11_tcells.py** - T cell subtype analysis
11. ✅ **step12_cell_communication.py** - Cell communication analysis

### Migration Strategy

1. **Sequential Migration**: Convert R scripts to Python in step-by-step order, each step as an independent executable script
2. **Centralized Configuration**: All file paths and parameters centralized in `config.py`, ensuring output structure matches R workflow
3. **Data Persistence**: Save intermediate AnnData objects (.h5ad) and CSV summaries for validation and debugging
4. **Quality Assurance**: Add lightweight regression tests to verify cluster counts, marker gene counts, and score distributions

## Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Run complete pipeline
python scripts/step01_qc_harmony.py
python scripts/step03_celltype.py
# ... run other steps sequentially
```

## Output Files

- **Data files**: `*.h5ad` (AnnData format, compatible with scanpy)
- **Statistical tables**: `*.csv` (marker genes, cell counts, statistical results)
- **Visualizations**: `*.png/*.pdf` (plots and heatmaps)
- **Logs**: `timings/*.jsonl` (runtime records for each step)

## Technical Highlights

- **Full Compatibility**: Results consistent with original R implementation
- **High Performance**: Optimized computational efficiency using numpy/scipy
- **Modular**: Each analysis step independent for easy debugging and extension
- **Standardized**: Uses scanpy standard data formats and best practices
- **Reproducible**: Fixed random seeds ensure reproducible results
