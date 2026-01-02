# scDEDS: Single-Cell Differential Equation Dynamical System for Gene Regulatory Network Inference

<img src="https://img.shields.io/badge/R%3E%3D-4.4.0-blue?style=flat&logo=R" /> <img src="https://img.shields.io/badge/Platform-Linux%20%7C%20Windows%20(WSL2)-lightgrey" /> <img src="https://img.shields.io/badge/License-MIT-yellow" />

## Overview
scDEDS is an R package designed for predicting Gene Regulatory Networks (GRNs) by modeling discrete pseudotemporal evolutionary dynamics from scRNA-seq and scATAC-seq data (the cells in the two datasets are identical; or, after using other methods to establish a one-to-one correspondence, the name of these matched cells is unified). By modeling gene expression dynamics as a **Differential Equation Dynamical System (DEDS)**, scDEDS can predict cell-specific regulatory interactions between Transcription Factors (TFs) and their Target Genes (TGs) across different cell states and pseudotemporal trajectories: TF-peak, peak-TG, TF-TG. It is linked to paper *scDEDS: A Discrete Evolutionary Dynamical System for GRN Inference from Paired Single-Cell Multi-omics Data*.

## Key Features

### 🧬 Multi-Omics Integration Framework
scDEDS pioneers a sophisticated integration approach that simultaneously leverages paired single-cell RNA-seq and ATAC-seq data. Unlike traditional methods that analyze these modalities separately, our package establishes direct connections between chromatin accessibility and gene expression dynamics, enabling more accurate inference of regulatory relationships.

### ⏱️ Pseudotemporal Dynamics Modeling
The core innovation of scDEDS lies in its differential equation-based modeling of gene regulatory dynamics along pseudotime trajectories. By treating cellular differentiation as a continuous dynamical system, we capture the temporal evolution of regulatory interactions, revealing how transcription factor activities and target gene responses change throughout developmental processes.

### 🌿 Branch-Aware Regulatory Analysis
scDEDS uniquely handles complex branching patterns in cell differentiation trajectories. Our algorithm automatically identifies branch points and constructs branch-specific gene regulatory networks, allowing researchers to investigate how regulatory programs diverge during lineage specification and cellular fate decisions.

### 🎯 TF Binding Site Integration
The package incorporates comprehensive transcription factor binding site information from the JASPAR database, combining computational predictions with empirical data to establish more reliable connections between transcription factors and their potential target genes.

### 📊 Multi-Layer Regulatory Inference
scDEDS provides a hierarchical regulatory network inference capability, predicting relationships at three interconnected levels:
- **TF → Chromatin fragments (peaks)**: TF-peak association predictions are obtained based on the motif binding scores of each transcription factor and each chromatin fragment.
- **Chromatin fragments (peaks) → Target Gene**: Peak-TG association predictions are obtained based on a human-defined promoter region (a base interval that includes the transcription start site, where open chromatin peaks within this interval are considered to have a regulatory relationship with the target gene corresponding to this transcription start site, such as promoters or enhancers).
- **TF → Target Gene**: One TG corresponds to several peaks. The match scores of transcription factors (TFs) to these peaks are calculated, and the average is taken as the initial regulatory strength of the TF for this TG. This process is performed for each TF and TG to obtain the initial gene regulatory network (iGRN). Typical TF-TG relationship pairs are selected as the training set to train the parameters of the DEDS model. For each TF-TG relationship pair, each set of trained parameters is applied, and the result with the best fit is selected as the predicted regulatory strength for that TF-TG. This process is iterated through for each TF-TG relationship pair to obtain the predicted gene regulatory network (pGRN).

## Full Workflow

### Phase 1: Data Preprocessing and Feature Engineering

**Target Gene Identification**
The workflow begins by identifying potential target genes through chromatin accessibility profiling. Promoter and enhancer regions are annotated using genome reference databases, linking accessible chromatin regions to their putative target genes based on genomic proximity and regulatory domain predictions.

**Transcription Factor Characterization**
Transcription factors are systematically cataloged using the JASPAR database, filtering for factors expressed in the cell type of interest. 

**Multi-Omics Data Integration**
Expression matrices for target genes and activity matrices for transcription factors are constructed, creating a unified data structure that synchronizes information from both transcriptomic and epigenomic modalities.

### Phase 2: Pseudotemporal Reconstruction and Cellular Organization

**Trajectory Inference**
Using dimensionality reduction and trajectory inference algorithms (monocle), cells are ordered along pseudotemporal axes that represent developmental or differentiation processes. This ordering captures the continuous nature of cellular transitions rather than treating cell states as discrete entities.

**Branch Point Identification**
The algorithm automatically detects branching points in the pseudotemporal trajectory (monocle), identifying where cell populations diverge into distinct lineages. Each branch represents an alternative developmental path with potentially unique regulatory programs.

**Dynamic Cell Grouping**
Cells are intelligently grouped along pseudotime, balancing temporal resolution with statistical power. The grouping strategy adapts to cellular density variations, ensuring sufficient cells in each temporal bin for robust statistical analysis while maintaining temporal precision.

### Phase 3: Initial Regulatory Network Construction

**TF-TG Relationship Scoring**
Based on transcription factor binding site predictions and co-expression patterns, initial regulatory relationships are scored. This creates a preliminary network that serves as the foundation for subsequent dynamical modeling.

**Branch-Specific Network Initialization**
Distinct initial networks are constructed for each identified branch, acknowledging that regulatory relationships may be branch-specific or exhibit differential strengths across lineages.

### Phase 4: Dynamical System Modeling and Parameter Optimization

**Differential Equation Framework**
The core innovation of scDEDS involves formulating gene regulation as a system of differential equations that describe the evolution in transcription factor expressions, target gene activities, and target gene expressions over pseudotime.

**Parameter Estimation via Hybrid Optimization**
A sophisticated optimization pipeline combining genetic algorithms and gradient descent methods estimates parameters that best explain the observed expression dynamics. This hybrid approach ensures robust convergence while avoiding suboptimal local minima.

**Model Selection and Validation**
Multiple candidate models are evaluated based on their ability to recapitulate observed expression patterns, with the best-performing models selected for final network prediction.

### Phase 5: Comprehensive Regulatory Network Prediction

**Final GRN Assembly**
The optimized dynamical models are used to predict regulatory strengths for all potential TF-TG pairs, generating a comprehensive, quantitative gene regulatory network.

**Multi-Layer Network Construction**
Beyond TF-TG interactions, the package infers relationships between TFs and chromatin accessibility, as well as between accessibility and gene expression, providing a multi-layer regulatory map.

**Quality Control and Filtering**
Stringent quality thresholds are applied to ensure high-confidence predictions, filtering out weak or unreliable interactions while retaining biologically meaningful regulatory relationships.

This comprehensive workflow transforms raw single-cell multi-omics data into dynamic, quantitative models of gene regulation, providing unprecedented insights into the mechanistic underpinnings of cellular differentiation and function.

## Hardware Recommendations

The model training process is **extremely computationally intensive**.
*   **OS:** Linux is highly recommended. Windows users can use WSL2.
*   **CPU:** A high-core-count CPU (e.g., >= 40 cores) is strongly advised.
*   **RAM:** A big RAM is highly recommended (e.g. for data with 1,000 cells and 20,000 genes, >= 64 GB is recommended, and more is of course better).
*   **Time:** Training for one cell type may take **several hours to days**.

## R Package

*   R package is in https://github.com/hth-buaa/scDEDS2.

## Specific Guidelines (Code, Data, and Results Availability)

*   Specific Guidelines (Code, Data, and Results Availability) is in https://github.com/hth-buaa/scDEDS-code-data-and-result/tree/main.
*   See specific Guidelines in *article_experiment_R_code.txt* (it is also the code for experiment in the paper).
*   See the experiment results of the paper in folder *article experiment result* (including the data used in paper). 
*   See the code for figures and some mediate analysis results in folder *article figure and some analysis result file R code*.
