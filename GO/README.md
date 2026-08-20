# Gene Ontology (GO) Analysis

This directory contains the scripts used to perform **Gene Ontology (GO) enrichment analyses** for the gene sets evaluated throughout the study.

The same workflow was run with different input gene sets to generate the GO analyses reported in the manuscript.

## Workflow

The analysis consists of two scripts:

```text id="oxmr8k"
Submission script that takes input gene list
      │
      ▼
GO enrichment analysis (R)
```

### Scripts

**Submission script**
Submits the GO analysis and specifies the input gene set and associated parameters.

**R analysis script**
Performs GO enrichment analysis for the supplied gene set and generates the corresponding results.

## Usage

To apply the workflow to a different gene set, modify the input specified by the submission script and rerun the analysis. The same underlying R workflow was used for the different GO analyses presented in the manuscript.

File paths and scheduler directives should be updated as necessary for the user's computing environment.

## Notes

* Input gene sets are derived from the corresponding analyses described in the manuscript.
* The workflow was originally executed in a high-performance computing environment.
