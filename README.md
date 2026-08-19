# asAPA_AD
Allele-specific alternative polyadenylation analysis on AD and control postmortem brains.

This repository contains the primary analysis and figure-generation code associated with our study of **Allele-specific alternative polyadenylation links noncoding genetic variation to Alzheimer's disease risk**.

The study investigates genetic regulation of alternative polyadenylation across human brain regions and integrates asAPA events with functional analyses.

The repository is organized around the major analyses presented in the manuscript. Because the complete analysis involved a large number of preprocessing, quality-control, and intermediate scripts, this repository focuses on the **core scripts required to reproduce the principal analyses described in the study**.

## Repository structure

```text
.
├── ASARP/
├── GO/
├── Concordance/
├── RBP/
├── FXS/
├── aQTL_GWAS/
├── REDIT/
└── README.md
```

### `ASARP/`

Scripts used to identify allele-specific alternative polyadenylation events from RNA-seq data using the ASARP framework and to process ASARP results for downstream analyses.

This directory contains code related to:

* preparation and processing of ASARP inputs and outputs
* identification of significant asAPA events
* aggregation of asAPA results across brain regions
* generation of summary statistics used in downstream analyses

### `GO/`

Gene Ontology enrichment analyses used to characterize biological processes and cellular functions associated with genes of interest.

This directory includes code for enrichment analyses performed on the principal asAPA gene set, prioritized functional subset, and FXS-shortened subset reported in the manuscript.

### `Concordance/`

Scripts used to evaluate the consistency of allele-specific effects across brain regions and prioritize variants exhibiting concordant asAPA behavior.

These analyses are used to distinguish reproducible cross-region regulatory effects and to define the high-concordance variant sets used in subsequent functional analyses.

### `RBP/`

Analyses evaluating potential RNA-binding protein mechanisms underlying asAPA regulation.

This directory contains code related to:

* eCLIP enrichment analyses
* sequence motif enrichment
* allele-specific RBP binding predictions

### `FXS/`

Scripts used to investigate alternative polyadenylation patterns in fragile X syndrome (FXS) datasets.

These analyses were used to provide additional evidence for a potential relationship between FMRP activity and transcript 3′-end regulation.

### `aQTL_GWAS/`

Code integrating asAPA results with genetic association datasets.

Analyses include:

* within-cohort alternative polyadenylation QTL (aQTL) mapping
* comparison of asAPA variants with aQTL signals
* integration with external aQTL resources
* overlap with genome-wide association study (GWAS) loci

### `REDIT/`

Scripts used to evaluate associations between allele-specific polyadenylation and Alzheimer’s disease phenotypes using **REDIT Regression**.

These analyses model allelic counts while accounting for relevant biological and technical covariates and were used to identify asAPA events associated with Alzheimer’s disease status and neuropathological traits.

## Analysis overview

The major analysis workflow can be summarized as:

```text
RNA-seq data
    │
    ▼
ASARP
Identification of asAPA events
    │
    ├──► Gene Ontology enrichment
    │
    ├──► Concordance
    │        │
    │        └──► RBP / motif analyses
    │
    ├──► aQTL and GWAS integration
    │
    └──► REDIT Regression
             │
             └──► Alzheimer’s disease associations

Additional FXS/FMR1 analyses
    │
    └──► Evaluation of FMRP-associated APA regulation
```

## Data

The analyses in this repository use human brain RNA-seq and genotype datasets described in the manuscript.

Due to data-use restrictions and the size of the underlying sequencing datasets, primary sequencing and genotype data are **not distributed directly through this repository**. Accession numbers and instructions for accessing the relevant public or controlled-access datasets are provided in the manuscript's **Data Availability** section.

Intermediate files required by individual scripts are described within the corresponding analysis directories where applicable.

## Software

The analyses were performed primarily using:

* Python
* R
* Bash
* ASARP
* REDIT Regression
* tensorQTL

Additional R and Python dependencies are specified within individual scripts or analysis directories.

Analyses were performed in a high-performance computing environment. Several scripts therefore contain or were originally designed for batch scheduling and may require modification of file paths, computational resources, or scheduler directives when executed in a different environment.

## Reproducibility

Scripts are organized according to the major analyses reported in the manuscript rather than as a single end-to-end executable workflow.

To reproduce a particular analysis:

1. Navigate to the corresponding analysis directory.
2. Review the script headers and required input files.
3. Update local input/output paths where necessary.
4. Install the required software and package dependencies.
5. Execute scripts in the order described within that directory.

Where analyses depend on intermediate outputs from another component of the study, the relevant dependency is noted in the corresponding script or directory documentation.

## Citation

If you use code from this repository, please cite:

> Barney, R. M., Quinones-Valdez, G., King, A. J., Amoah, K., Wang, W., & Xiao, X. (2026). Allele-specific alternative polyadenylation links noncoding genetic variation to Alzheimer's disease risk. bioRxiv : the preprint server for biology, 2026.02.13.705798. https://doi.org/10.64898/2026.02.13.705798

Citation information will be updated following publication.

## Code availability

The code in this repository accompanies the analyses described in the associated manuscript and is provided to facilitate transparency and reproducibility of the reported results.

Questions regarding the analyses or implementation can be submitted through the repository's GitHub issue tracker.
