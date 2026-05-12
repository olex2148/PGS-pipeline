# PGS-pipeline
Pipeline using [gwf workflow](https://gwf.app/) and [LDpred2-auto](https://privefl.github.io/bigsnpr/articles/LDpred2.html) 

This repository contains a scalable and reproducible pipeline for harmonising GWAS summary statistics and estimating key parameters of genetic architecture using LDpred2-auto. The pipeline was developed for the study:

Hansen et al. Characterising the genetic architecture of thousands of complex traits from public GWAS summary statistics.

Overview

The pipeline enables systematic estimation of:

SNP-based heritability (h²)
Polygenicity (p)
Frequency-dependent effect size parameter (α)
Inferred predictive performance (R²)

It operates on publicly available GWAS summary statistics (e.g. from the GWAS Catalog) and applies standardised quality control, harmonisation, and model fitting within a unified framework.

Requirements: 

Linux 
Conda environment with the following
GWF workflow (https://gwf.app/) ≥  v2.1 
Python ≥ v3.8
R ≥ v4.2

R packages: bigsnpr, bigstatsr, MungeSumstats, data.table, rjson, openxlsx, testit, dplyr, stringr, ggplot2

Relevant folder structure: 

.gwf: 
	Here the logfiles output by gwf can be found, useful for debugging
code: 
	Contains scripts that munge and generates PGS model. This is also where "n" can be added if it is missing in the sumstats. 
	Also contains the "functions" folder that rename sumstat columns, calculates missing variables depending on what is present in the sumstats, extracts n, converts OR to beta, etc 
data: 
	data contains HapMap3+, the LD reference, and a folder called "sumstats" where all sumstats should be placed. 
	In the original pipeline it contains paths.json for path specification to genotype data, HapMap3+, phenotype data, LD blocks, iPSYCH meta data, relatives and pcs (used for making a homogenous population when calculating adjusted R2). In this standalone pipeline run these paths are modified to be generalised for future users and this test. The json also is where one should specify ones work directory "work_dir" to the desired subfolder in "sumstats" according to the sumstats that are desired to be run. 
root: 
	Where the workflow.py is (i.e. ../data)
results: 
	Here the outputs will automatically be placed
steps: 
	The munged sumstats will be placed here

Instructions/Commands: 

When the required packages and software is installed and activated in ones conda environment and the user is situated in the root folder (i.e. where the workflow.py is) one can give the "gwf run" command in the linux terminal. This will submit the jobs that have not yet run in the specified working directory. 
Another relevant gwf command is "gwf status" which will let the user know if the munging or model is complete or still running, or if it has failed. If it has failed go to the .gwf folder and inspect the logs. Most of the time "n" is missing and needs to be added manually in add_n from the code folder. 
 

Outputs:

Logfile in .gwf folder
Variant QC plot
1st kept (converged) chain
LDpred2-auto PGSs
Lassosum2-auto PGSs
Meta-data sheet relevantly including model parameters h2, inf-r2, polygenicity, selection coef alpha, LDSC h2, variant count before and after qc, number of kept chains, information on the study (year, title, authors, accession id), n_eff, cases, controls, shrinkage, trait, misc

The pipeline is designed for large-scale application (>1,000 traits) and ensures consistent parameter estimation across heterogeneous GWAS datasets.
All filtering thresholds and QC criteria are implemented programmatically and documented in the manuscript.
Users are responsible for ensuring compatibility between GWAS ancestry and the LD reference panel.
Citation

If you use this pipeline, please cite the associated manuscript (details to be updated upon publication).
