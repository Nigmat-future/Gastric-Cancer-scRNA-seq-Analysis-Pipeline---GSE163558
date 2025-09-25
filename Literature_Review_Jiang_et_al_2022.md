# Literature Review: Revealing the transcriptional heterogeneity of organ-specific metastasis in human gastric cancer using single-cell RNA Sequencing

## 1. Core Information

- **Title**: Revealing the transcriptional heterogeneity of organ-specific metastasis in human gastric cancer using single-cell RNA Sequencing
- **Journal**: *Clinical and Translational Medicine* (2022)
- **DOI**: `10.1002/ctm2.730`
- **Research Fields**: Gastric Cancer, Tumor Metastasis, Single-Cell RNA Sequencing, Tumor Microenvironment

---

## 2. Background and Objective

Gastric cancer is one of the most common malignancies worldwide, with a high mortality rate primarily due to tumor metastasis. The microenvironments of different organs (e.g., liver, peritoneum, lymph nodes) vary significantly, but the molecular mechanisms by which tumor cells adapt to these specific environments to achieve successful metastasis are not well understood.

This study aims to use **single-cell RNA sequencing (scRNA-seq)** technology to dissect the **transcriptional heterogeneity** of gastric cancer across the **Primary Tumor (PT)**, **Lymph Node metastasis (LN)**, and **Peritoneal metastasis (P)**, in order to uncover the cellular and molecular mechanisms of organ-specific metastasis.

---

## 3. Methods

- **Sample Source**: 10 tissue samples were collected from 6 gastric cancer patients, including primary tumors (PT), lymph node metastases (LN), and peritoneal metastases (P).
- **Technology Platform**: **10x Genomics Chromium** platform was used for single-cell sequencing.
- **Data Analysis Pipeline**:
    - **Quality Control (QC)**: Low-quality cells were filtered out (genes < 200 or > 5000, mitochondrial gene ratio > 20%).
    - **Data Integration and Dimensionality Reduction**: The `Seurat` package was used for data normalization, batch correction, PCA, and t-SNE/UMAP visualization.
    - **Cell Type Annotation**: Cells were annotated based on known marker genes, primarily classified into epithelial cells, immune cells, and stromal cells.
    - **Malignant Cell Identification**: `inferCNV` was used to infer copy number variations in epithelial cells to identify malignant tumor cells.
    - **Differential Expression Analysis**: Gene expression differences were compared among malignant cells from different metastatic sites (PT, LN, P).
    - **Functional Enrichment Analysis**: GO and KEGG enrichment analyses were performed on differentially expressed genes.
    - **Cell-Cell Communication Analysis**: `CellChat` was used to analyze interactions between different cell subpopulations within the tumor microenvironment.

---

## 4. Key Findings

### Finding 1: Gastric cancer metastases exhibit high cellular heterogeneity
- The study successfully identified 8 major cell types, including epithelial cells, T cells, B cells, macrophages, and fibroblasts.
- The cellular composition varied significantly across different metastatic sites. For example, peritoneal metastases had a higher proportion of fibroblasts, while lymph node metastases were enriched with immune cells.

### Finding 2: Malignant cells display unique transcriptional signatures in different metastatic organs
- **Peritoneal Metastasis (P)**:
    - Malignant cells showed high expression of genes related to **Epithelial-Mesenchymal Transition (EMT)**, **Extracellular Matrix (ECM) remodeling**, and **angiogenesis** (e.g., `TGFBI`, `COL1A1`).
    - This suggests that tumor cells adapt to the unique peritoneal microenvironment by activating these pathways.
- **Lymph Node Metastasis (LN)**:
    - Malignant cells upregulated pathways associated with **immune evasion**, such as the downregulation of **MHC class I molecules** and expression of **PD-L1** (`CD274`).
    - This helps tumor cells survive in the immune-rich environment of the lymph nodes.
- **Primary Tumor (PT)**:
    - Malignant cells in the primary tumor exhibited stronger **proliferative properties**, with higher expression of cell-cycle-related genes (e.g., `MKI67`).

### Finding 3: The tumor microenvironment plays a crucial role in organ-specific metastasis
- **Peritoneal Metastasis Microenvironment**:
    - A specific subtype of **Cancer-Associated Fibroblasts (CAFs)** was identified. These CAFs interact with tumor cells by secreting cytokines (e.g., `TGF-β`), collectively promoting peritoneal seeding.
- **Lymph Node Metastasis Microenvironment**:
    - Strong interactions were observed between tumor cells and **Exhausted T cells**, indicating the formation of an immunosuppressive microenvironment.

### Finding 4: Identification of potential therapeutic targets
- The study found that the chemokine receptor **CCR1** was specifically overexpressed in malignant cells of peritoneal metastases. In vitro experiments confirmed that inhibiting `CCR1` could effectively reduce the invasive ability of tumor cells, suggesting `CCR1` as a potential therapeutic target for treating peritoneal metastasis in gastric cancer.

---

## 5. Conclusion and Significance

This study provides the first systematic single-cell transcriptomic atlas of gastric cancer metastasis across different organs, revealing how **tumor cells** and the **microenvironment** synergistically drive organ-specific metastasis.

**Key Conclusions**:
1.  Gastric cancer metastasis is not a random process but involves **programmatic transcriptional reprogramming** of tumor cells to adapt to specific organ microenvironments.
2.  **Peritoneal metastasis** is primarily associated with EMT and matrix remodeling, whereas **lymph node metastasis** is closely linked to immune evasion mechanisms.
3.  The identification of organ-specific genes (like `CCR1`) offers new avenues for developing targeted therapies for specific types of metastasis.

## 6. Implications for the Current Project (GSE163558)

This paper provides an excellent reference framework and theoretical support for your ongoing gastric cancer project:

- **Validate Your Findings**: You can compare whether the analysis results from your `GSE163558` dataset (which includes more metastatic types like liver and ovary) are consistent with the findings of this paper.
- **Adopt Analysis Methods**: The methods used in this paper, such as `inferCNV` for malignant cell identification and `CellChat` for communication analysis, can be directly applied to your project.
- **Extend Your Analysis**:
    - You can focus on the specific genes and pathways in the **liver (Li)** and **ovary (O)** metastasis samples from your data.
    - Compare the **immune microenvironment** differences across various metastatic sites (LN, P, Li, O), especially the exhaustion status of T cell subpopulations.
    - Identify "common" metastatic genes that are upregulated or downregulated across all metastatic sites, as well as "organ-specific" genes that change only in particular organs.
