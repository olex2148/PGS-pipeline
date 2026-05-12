# PGS-pipeline
Pipeline using [gwf workflow](https://gwf.app/) and [LDpred2-auto](https://privefl.github.io/bigsnpr/articles/LDpred2.html) 


This repository contains a scalable and reproducible pipeline for harmonising GWAS summary statistics and estimating key parameters of genetic architecture using LDpred2-auto. The pipeline was developed for the study:



Hansen et al. Characterising the genetic architecture of thousands of complex traits from public GWAS summary statistics.

# IMPORTANT <br>

Download [map_hm3_plus.rds](https://figshare.com/articles/dataset/LD_reference_for_HapMap3_/21305061) and place it in the data/hapmap3+ folder. Then download [ldref_with_blocks.zip](https://figshare.com/articles/dataset/European_LD_reference_with_blocks_/19213299) and unzip in the data/ld_blocks folder.

## Overview

The pipeline enables systematic estimation of:

- SNP-based heritability (h²)
- Polygenicity (p)
- Frequency-dependent effect size parameter (α)
- Inferred predictive performance (R²)

It operates on publicly available GWAS summary statistics (e.g. from the GWAS Catalog) and applies standardised quality control, harmonisation, and model fitting within a unified framework.

## Requirements: 

- Linux 
- Conda environment with the following <br>
- GWF workflow (https://gwf.app/) ≥  v2.1 <br>
- Python ≥ v3.8<br>
- R ≥ v4.2<br>

## R packages: <br>

**bigsnpr, bigstatsr, MungeSumstats, data.table, rjson, openxlsx, testit, dplyr, stringr, ggplot2** <br>

## Relevant folder structure: <br>

.gwf: <br>
	Here the logfiles output by gwf can be found, useful for debugging<br>
code: <br>
	Contains scripts that munge and generate the PGS model. This is also where "n" can be added if it is missing in the sumstats. <br>
	Also contains the "functions" folder that rename sumstat columns, calculates missing variables depending on what is present in the sumstats, extracts n, converts OR to beta, etc <br>
data: <br>
	data contains HapMap3+, the LD reference, and a folder called "sumstats" where all sumstats should be placed. <br>
	In the original pipeline it contains paths.json for path specification to genotype data, HapMap3+, phenotype data, LD blocks, iPSYCH meta data, relatives and pcs (used for making a homogenous population when calculating adjusted R2). In this standalone pipeline run these paths are modified to be generalised for future users and this test.  <br>
root: <br>
	Where the workflow.py is (i.e. ../data)<br>
results: <br>
	Here the outputs will automatically be placed<br>
steps: <br>
	The munged sumstats will be placed here<br>

## Instructions/Commands: <br>

Firstly download and place the LD block and Hapmap3+ files as described above. Next edit the json in "data" and specify the path to your work directory "work_dir" e.g. ~/C:/Users/USER/X/Desktop/PGS_workflow/. 

When the required packages and software are installed and activated in ones conda environment and the user is situated in the root folder (i.e. where the workflow.py is) one can give the "gwf run" command in the linux terminal. This will submit the jobs that have not yet run in the specified working directory. <br>
Another relevant gwf command is "gwf status" which will let the user know if the munging or model is complete or still running, or if it has failed. If it has failed go to the .gwf folder and inspect the logs. Most of the time "n" is missing and needs to be added manually in add_n from the code folder. <br>
 

## Outputs:<br>

Logfile in .gwf folder<br>
Variant QC plot<br>
1st kept (converged) chain<br>
LDpred2-auto PGSs<br>
Lassosum2-auto PGSs<br>
Meta-data sheet relevantly including model parameters h2, inf-r2, polygenicity, selection coef alpha, LDSC h2, variant count before and after qc, number of kept chains, information on the study (year, title, authors, accession id), n_eff, cases, controls, shrinkage, trait, misc<br>

All filtering thresholds and QC criteria are implemented programmatically and documented in the manuscript and should be adapted to fit the data that is being worked on.<br>
Users are responsible for ensuring compatibility between GWAS ancestry and the LD reference panel.<br>
Installation time: N/A
Special hardware: None
Runtime: Depends on sumstat size, test set takes ~10 min 
## Citation

If you use this pipeline, please cite the associated manuscript (details to be updated upon publication).
