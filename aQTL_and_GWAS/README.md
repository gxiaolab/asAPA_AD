# aQTL and GWAS Analysis

This directory contains the scripts used to integrate allele-specific alternative polyadenylation (asAPA) results with **alternative polyadenylation quantitative trait loci (aQTLs)** and **genome-wide association study (GWAS)** signals.

The workflow consists of two analysis modules: within-cohort aQTL mapping using **tensorQTL** and linkage disequilibrium (LD)-based comparison of asAPA variants with GWAS-associated variants.

## Workflow

Scripts are numbered according to their execution order and are divided into two analysis modules:

```text id="trwn7z"
asAPA and genotype data
          │
          ├── 1–8. aQTL analysis
          │       │
          │       ├── Prepare brain-region APA matrices
          │       ├── Prepare tensorQTL inputs
          │       └── Run tensorQTL and process results
          │
          └── 9–10. GWAS analysis
                  │
                  ├── Prepare variants for LD analysis
                  └── Identify LD with GWAS variants
```

Required helper scripts for individual analysis steps are included in this directory.

## Analysis modules

### Scripts 1–8 — aQTL analysis

These scripts perform within-cohort aQTL mapping for the brain regions analyzed in the study.

**Scripts 1–2 — APA matrix preparation**
Generate and organize the brain-region-specific alternative polyadenylation matrices required for QTL mapping.

**Scripts 3–4 — tensorQTL input preparation**
Prepare APA phenotypes, genotype information, and associated inputs in the formats required by tensorQTL.

**Scripts 5–8 — tensorQTL analysis and output processing**
Run tensorQTL to identify genetic variants associated with alternative polyadenylation and process the resulting aQTL statistics for downstream comparison with asAPA variants.

### Scripts 9–10 — GWAS linkage disequilibrium analysis

These scripts evaluate relationships between asAPA variants and disease-associated variants identified through GWAS.

**Script 9 — LD analysis preparation**
Prepares asAPA and GWAS variants for linkage disequilibrium analysis.

**Script 10 — GWAS LD analysis**
Evaluates linkage disequilibrium between asAPA variants and GWAS-associated variants and organizes the resulting disease-associated loci for downstream analysis.

## External tools and data sources

### tensorQTL

Within-cohort aQTL mapping was performed using **tensorQTL**, a GPU-enabled QTL mapping framework.

**tensorQTL:** (https://github.com/broadinstitute/tensorqtl)

The software version, QTL mapping parameters, covariates, and statistical criteria used for the study are described in the accompanying manuscript and relevant scripts.

### GWAS data

GWAS variants and summary statistics were obtained from previously published studies described in the accompanying manuscript. This data can be accessed on the GWAS catalog.

**GWAS Catalog:** (https://www.ebi.ac.uk/gwas/)

## Running the analyses

Scripts within each module should be executed according to their numerical prefixes. 

The aQTL and GWAS modules have distinct external input requirements. File paths, genotype resources, covariates, and computational parameters should therefore be reviewed and modified as necessary for the user's computing environment.

## Notes

* aQTL mapping was performed separately for the brain regions analyzed in the study.
* tensorQTL results represent within-cohort associations between genetic variation and alternative polyadenylation phenotypes.
* GWAS analyses evaluate genetic relationships through linkage disequilibrium and do not imply that an asAPA tag SNP is itself the causal disease-associated variant.
* The workflow was originally executed in a high-performance computing environment and may require modification for other computational environments.
