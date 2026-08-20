# Fragile X Syndrome (FXS) Analysis

This directory contains the analyses used to investigate alternative polyadenylation patterns in **Fragile X syndrome (FXS)** and their relationship to FMRP-associated sequence motifs.

The workflow compares transcript 3′-end usage between FXS and control samples to identify relative transcript lengthening or shortening and subsequently evaluates whether these changes are associated with differences in FMRP motif density.

## Workflow

The analysis consists of **five sequential scripts**, numbered according to their required execution order:

```text id="c7l2yf"
PolyA_DB annotations
        │
        ▼
1–2. Prepare polyadenylation site annotations
        │
        ▼
3. Compare FXS and control samples
   Identify transcript lengthening / shortening
        │
        ▼
4. Prepare sequences for motif analysis
        │
        ▼
5. Compare FMRP motif density between
   lengthened and shortened transcripts
```

### Analysis stages

**Scripts 1–2 — Polyadenylation site preparation**
Process and organize polyadenylation site annotations obtained from **PolyA_DB** for use with the FXS datasets.

**Script 3 — Transcript length analysis**
Compares polyadenylation site usage between FXS and control samples within each cohort to determine whether transcripts exhibit relative 3′ UTR lengthening or shortening in FXS compared to controls.

**Script 4 — Motif analysis preparation**
Prepares the sequences and transcript groups required for downstream FMRP motif analysis.

**Script 5 — FMRP motif density analysis**
Compares the density of FMRP-associated sequence motifs between transcripts exhibiting lengthening and shortening to evaluate whether FXS-associated changes in polyadenylation are associated with potential FMRP regulatory sites.

## Data sources

Polyadenylation site annotations were obtained from **PolyA_DB (https://exon.apps.wistar.org/polya_db/v3/)**.

The FXS RNA-seq datasets analyzed here were obtained from previously published cohorts. Dataset accession information and original study references are provided in the accompanying manuscript's **Data Availability** section.

## Running the analysis

Scripts should be executed in numerical order (`1`–`5`).

File paths should be updated to reflect the user's local data organization and computing environment.

## Notes

* PolyA_DB annotations are not redistributed here and should be obtained from the original resource.
* See the accompanying manuscript for complete descriptions of the FXS cohorts, statistical analyses, and interpretation of these results.
