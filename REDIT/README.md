# REDIT Regression Analysis

This directory contains the scripts used to test associations between **allele-specific alternative polyadenylation (asAPA)** and Alzheimer’s disease phenotypes using **REDIT Regression**.

The workflow prepares allele-specific count data and associated sample metadata for regression analysis, runs a study-specific adaptation of REDIT Regression, and organizes the resulting statistical output for downstream analysis.

## Workflow

The analysis consists of two primary scripts and an adapted REDIT Regression R script:

```text id="r8q1md"
1. Prepare REDIT Regression input (obtain allele-specific counts)
                 │
                 ▼
2. Run REDIT Regression and organize results
```

### Script 1 — Input preparation

Prepares the allele-specific count data, phenotype information, and covariates required for REDIT Regression.

### Script 2 — Regression and output processing

Submits the REDIT Regression analysis using the prepared inputs and organizes the resulting statistical output for downstream analyses and visualization.

### Adapted REDIT Regression script

A study-specific adaptation of the REDIT Regression R script is included in this directory. This version was adapted for the allele-specific polyadenylation analyses and phenotypes evaluated in this study.

## REDIT Regression

REDIT Regression was originally developed by the Xiao Lab. The original implementation and documentation are publicly available at:

**(https://github.com/gxiaolab/REDITs)**

The adapted implementation included here reflects the version used for the analyses reported in the accompanying manuscript. Readers interested in the original framework and complete documentation should refer to the repository above.

## Running the analysis

File paths, phenotype variables, covariates, and scheduler directives should be reviewed and modified as necessary for the user's dataset and computing environment.

## Notes

* The included REDIT Regression R script has been adapted for the analyses performed in this study and should not be considered the canonical distribution of REDIT Regression.
* Primary sequencing, genotype, and phenotype data are not distributed with this repository. See the manuscript's **Data Availability** section for information on accessing the datasets used in the study.
* The workflow was originally executed in a high-performance computing environment and may require modification for other computational environments.
