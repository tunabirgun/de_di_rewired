# Differential Interactome & Rewiring Workflow

This repository contains the computational pipeline used to analyze differential protein-protein interactions (PPI) and identify functional rewiring events. The workflow integrates RNA-Seq expression data with static PPI networks to infer condition-specific interactions.

## References

**Original Method**: The differential interactome algorithm (Q-value scoring based on binarized expression) is an implementation of the method described in:

> Gulfidan, G., Turanli, B., Beklen, H. *et al.* (2020). "**Pan-cancer mapping of differential protein-protein interactions**". *Scientific Reports*. 10, 3272.
> [https://doi.org/10.1038/s41598-020-60127-x](https://www.nature.com/articles/s41598-020-60127-x)

Please cite this original paper when using the differential interactome component of this pipeline.

## Pipeline Overview

The analysis follows a six-stage protocol:

0.  **Quantification (`00`)**: RNA-Seq analysis pipeline (SRA -> FASTQ -> BAM -> Counts).
1.  **Differential Expression (`01`)**: Transcriptome-level analysis (DESeq2) to identify DEGs.
2.  **Preprocessing (`02`)**: Standardization of gene IDs and orthology mapping.
3.  **Differential Analysis (`03`)**: The core algorithm detecting interaction changes based on expression states.
4.  **Rewiring Extraction (`04`)**: Filtering interactions to isolate "Rewired" proteins.
5.  **Enrichment (`05`)**: Functional characterization of rewired candidates.
6.  **Visualization (`06`)**: Generation of the figures.

## Script Usage

### 0. `00_rnaseq_quantification.sh`
A shell script template for processing raw SRA data.
*   **Tools**: `sra-tools`, `fastp`, `SortMeRNA`, `STAR`, `featureCounts`.

### 1. `01_differential_expression.R`
**Purpose**: Standard differential expression analysis.
*   **Method**: `DESeq2` (Normalization, Wald Test, LFC Shrinkage).
*   **Outputs**:
    *   `_Upregulated.xlsx` / `_Downregulated.xlsx` (Used as input for Step 04).
    *   **Figures**: PCA, Volcano Plots, Heatmaps (Top 50 DEGs).

### 2. `02_preprocess_ids.R`
Standardizes input matrices to a common Gene ID format required by the downstream R scripts.

### 3. `03_differential_analysis.R`
**Core Method**: Implements the Differential Interactome algorithm.
*   **Source**: Adjusted implementation of *Gulfidan et al., Sci Rep (2020)*.
*   **Binarization**: Converts continuous expression counts to discrete states (1, 0, -1) based on local mean and standard deviation.
*   **Q-Value**: Calculates the probability of an interaction occurring in the Treatment vs Control group.
*   **Classification**: Classifies interactions as "Activated", "Repressed", or "Unchanged".

### 4. `04_rewired_analysis.R`
Separates "Abundance-Driven" changes from "True Rewiring".
*   **Logic**: Filters for proteins with high differential degree (top 10%) but stable expression (Low LogFC).
*   **Output**: A curated list of rewired candidates.

### 5. `05_enrichment.R`
Performs Gene Ontology (GO) enrichment analysis on the rewired gene lists using `clusterProfiler`.

### 6. `06_visualization.R`
Generates figures:
*   **UpSet Plots**: To visualize intersections of rewired proteins across multiple conditions.
*   **Venn Diagrams**: To compare specific subsets (e.g., Gained vs Lost interactions).

---
*Note: This workflow requires R 4.0+ and standard Bioconductor packages (`dplyr`, `AnnotationHub`, `clusterProfiler`, `UpSetR`).*
