# BCR Repertoire Analysis

This repository contains the code and analysis pipeline for my Master's thesis at Aarhus Universitet. The project investigates features of the B-cell receptor (BCR) immune repertoire, including diversity, clonality, gene usage, and clinical correlations, in **renal** and **bladder** cancer patients, compared against healthy controls. The goal is to better understand the BCR’s role in anti-cancer immunity and its potential clinical relevance.


## Data

The project uses preprocessed BCR repertoire data stored as RDS objects (`renal_data`, `bladder_data`) with MiXCR-derived fields:
`mixcr_df_IGH`, `mixcr_df_IGK`, `mixcr_df_IGL`

Clinical metadata is integrated within these RDS files.

> **Note**: Raw data is not included due to patient confidentiality. 

## Repertoire Analyses

### 1. Clonotype and Expansion
- **File**: `clonotype_analysis.R`
- Outputs boxplots of unique clonotypes, read depth, and clonal expansion profiles.

### 2. Diversity Metrics
- **File**: `diversity_analysis.R`
- Computes and visualizes:
  - Shannon & Normalized Shannon
  - Gini Index
  - Clinical associations with normalized diversity

### 3. Survival Analysis
- **File**: `survival_analysis.R`
- Kaplan-Meier survival stratified by Shannon diversity
- Cox proportional hazards models with age adjustment

### 4. Healthy vs Cancer Comparison
- **File**: `cancerVShealthy_analysis.R`
- Compares read depth, diversity, and clonal expansion between healthy, renal, and bladder donors

### 5. UMAP of Gene Usage
- **File**: `umap.R`
- Visualizes gene usage patterns using UMAP (colored by age, sex, health status)

### 6. Gene Usage Correlation with Clinical Features
- **File**: `gene_usage_correlation_metadata.R`
- Computes and plots Spearman correlations (heatmaps)

## Requirements

- R ≥ 4.2
- R packages:
  - `tidyverse`, `vegan`, `ggpubr`, `survival`, `survminer`, `ComplexHeatmap`, `umap`, `pheatmap`, etc.

---

**Thesis Title**: Investigating the B-Cell Receptor in bladder and kidney cancer to understand its role in anti-cancer immunity and its potential implications for cancer treatment and outcomes  
**Author**: Liliane Zoe Bader  
**Submitted**: June 15, 2025  
**Program**: Master’s in Bioinformatics  
**Institution**: Aarhus University, Denmark  
**Department**: Bioinformatics Research Center (BiRC)  
**Thesis Host Lab**: Cancer Evolution & Immunology Group, Department of Molecular Medicine (MoMa), Aarhus University Hospital  
**Supervisor**: Prof. Nicolai Juul Birkbak  

