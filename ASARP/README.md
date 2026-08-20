# ASARP Analysis

This directory contains the core scripts used to identify and process allele-specific alternative polyadenylation (asAPA) events from RNA-seq data using **ASARP**.

## Workflow

The analysis consists of two main stages:

1. **RNA-seq processing and ASARP analysis**
   Performs the primary processing workflow, including:

   * alignment with STAR
   * marking and removing duplicates
   * finding RNA-DNA differences
   * asAPA detection with ASARP

2. **ASARP output processing**
   Organizes and summarizes the raw ASARP output for downstream analyses. Python and R helper scripts used during this processing step are included in this directory.

## Scripts

```text
ASARP/
├── 1_asarp.sh
├── 2_asarp_organize.sh
├── organize1.py
└── organize2.py
```

## ASARP and RDD

ASARP and RDD were developed and maintained by the Xiao Lab. The complete implementations and documentation for these tools are maintained separately in the lab's private GitHub repositories and therefore are not distributed as part of this repository.

The scripts provided here contain the project-specific workflow used to apply these methods to the datasets analyzed in this study and the subsequent processing of their outputs.

## Notes

* File paths should be updated to match the user's local computing environment.
* The workflow was originally executed in a high-performance computing environment and may require modification of scheduler directives and computational resources on other systems.
* Primary sequencing data are not included in this repository. See the manuscript's **Data Availability** section for information on accessing the datasets used in the study.
