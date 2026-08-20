# Concordance Analysis

This directory contains the pipeline used to identify **putative functional SNPs** based on concordance of allele-specific alternative polyadenylation (asAPA) signals across individuals.

The workflow evaluates allele-specific read distributions surrounding candidate variants, quantifies concordance between allelic signals, and uses probabilistic modeling and multiple-testing correction to prioritize variants with evidence of functional effects on alternative polyadenylation.

## Workflow

The analysis is implemented as a series of **sequential Bash scripts**, numbered according to their required execution order. Python helper scripts called by individual steps are included in the same directory.

At a high level, the pipeline performs the following:

```text id="l2yxak"
Candidate asAPA regions
        │
        ▼
Extract alternative 3′ UTR regions
        │
        ▼
RDD processing
        │
        ▼
Candidate organization and VCF filtering
        │
        ▼
Calculate concordance scores
        │
        ▼
Gaussian mixture model
Estimate concordance probabilities
        │
        ▼
Classify putative functional SNPs
(v1 / v2 / v2b)
        │
        ▼
FDR correction
        │
        ▼
Summarize prioritized functional SNPs
```

### Major analysis stages

The numbered scripts collectively perform:

1. extraction and preparation of alternative 3′ UTR regions containing candidate asAPA variants
2. RDD processing of the selected regions
3. organization and filtering of candidate variants using genotype information from VCF files
4. calculation of allele-specific concordance scores
5. estimation of concordance probabilities using a Gaussian mixture model
6. classification of candidate SNPs into putative functional categories (`v1`, `v2`, and `v2b`)
7. false discovery rate (FDR) correction
8. generation and summarization of the final prioritized SNP set

## Running the pipeline

Bash scripts are numbered `1` through `21` and should be executed **in numerical order**.

Python scripts in this directory are helper scripts invoked by the corresponding Bash workflow steps and generally should not be run independently.

File paths, computational resources, and scheduler directives should be modified as necessary for the user's computing environment.

## Functional SNP classifications

The pipeline assigns prioritized variants to the `v1`, `v2`, and `v2b` categories according to the concordance-based criteria implemented by the original method.

These classifications represent distinct patterns of allele-specific evidence used to identify variants with putative functional effects. See the original concordance pipeline and associated methodology for the complete definitions and statistical framework underlying these classifications.

## Original implementation

The concordance-based functional SNP prioritization pipeline was originally developed by the Xiao Lab and is publicly available in the lab's GitHub repository:

**(https://github.com/gxiaolab/cGMAS)**

The scripts provided in this directory represent the study-specific implementation of this pipeline used for the analyses reported in the accompanying manuscript.

## Notes

* Primary sequencing and genotype data are not distributed with this repository. See the manuscript's **Data Availability** section for information on accessing the relevant datasets.
* RDD is a dependency of portions of this workflow; its complete implementation and documentation are maintained separately by the Xiao Lab.
* The pipeline was originally executed in a high-performance computing environment and may require modification for other computational environments.
