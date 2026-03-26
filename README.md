# PGS-pipeline
Pipeline using [gwf workflow](https://gwf.app/) and [LDpred2-auto](https://privefl.github.io/bigsnpr/articles/LDpred2.html) 

This repository contains a scalable and reproducible pipeline for harmonising GWAS summary statistics and estimating key parameters of genetic architecture using LDpred2-auto. The pipeline was developed for the study:

Hansen et al. Characterising the genetic architecture of thousands of complex traits from public GWAS summary statistics (under review).

Overview

The pipeline enables systematic estimation of:

SNP-based heritability (h²)
Polygenicity (p)
Frequency-dependent effect size parameter (α)
Inferred predictive performance (R²)

It operates on publicly available GWAS summary statistics (e.g. from the GWAS Catalog) and applies standardised quality control, harmonisation, and model fitting within a unified framework.

Requirements
R (≥ 4.2 recommended)
Key packages: bigsnpr, bigstatsr, MungeSumstats, data.table
Access to an LD reference panel (e.g. HapMap3+)
Optional: HPC environment for parallel execution across traits
Running the pipeline
Prepare GWAS summary statistics
Input files should contain effect sizes and standard errors. The pipeline includes harmonisation (via MungeSumstats) and performs variant-level QC automatically.

Specify input paths
Paths to:

GWAS summary statistics
LD reference panel
(Optional) target genotype data for PGS evaluation

Note: All file paths have been removed from the repository and must be specified by the user prior to execution.

Run main script
The pipeline can be executed trait-wise or in parallel across multiple GWAS datasets. 

Outputs
Results are saved as .rds files and include:
Genetic architecture estimates (h², p, α)
Model diagnostics and convergence metrics
QC summaries
(Optional) polygenic scores
Notes on reproducibility
The pipeline is designed for large-scale application (>1,000 traits) and ensures consistent parameter estimation across heterogeneous GWAS datasets.
All filtering thresholds and QC criteria are implemented programmatically and documented in the manuscript.
Users are responsible for ensuring compatibility between GWAS ancestry and the LD reference panel.
Citation

If you use this pipeline, please cite the associated manuscript (details to be updated upon publication).
