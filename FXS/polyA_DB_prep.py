#!/usr/bin/python
import sys
import re
import pandas as pd
import numpy as np
from itertools import groupby

## Get positions

file_path = sys.argv[1]
pas_path = sys.argv[2]
endt_path = sys.argv[3]

# Define column names
column_names = ["chrom", "PAS_pos", "strand", "gene_id", "gene_number", "gene_name", "PAS_ID",
                "transcript_id", "transcript_start", "transcript_end", "strand_dup", "utr_length", "region"]

# Load data
df = pd.read_csv(file_path, sep="\t", header=None, names=column_names)

# Keep relevant columns and convert types
df = df[["chrom", "PAS_pos", "transcript_start", "transcript_end", "strand"]]
df["PAS_pos"] = pd.to_numeric(df["PAS_pos"], errors="coerce")
df["transcript_start"] = pd.to_numeric(df["transcript_start"].astype(str).str.replace('.0', '', regex=False), errors="coerce")
df["transcript_end"] = pd.to_numeric(df["transcript_end"].astype(str).str.replace('.0', '', regex=False), errors="coerce")

# Drop any rows with missing values
df.dropna(inplace=True)

# Drop rows where PAS and transcript site are too close (500bp), strand-aware
def is_valid_distance(row):
    if row["strand"] == "+":
        return abs(row["transcript_end"] - row["PAS_pos"]) >= 400
    else:
        return abs(row["transcript_start"] - row["PAS_pos"]) >= 400

# Drop rows where PAS and transcript end are < 500bp apart (absolute distance)
df = df[df.apply(is_valid_distance, axis=1)]

# Create strand-aware PAS windows
df["PAS_window"] = np.where(df["strand"] == "+",
                            df["PAS_pos"] - 400,
                            df["PAS_pos"] + 400)

# Create strand-aware transcript windows
df["transcript_window"] = np.where(df["strand"] == "+",
                                   df["transcript_end"] - 400,
                                   df["transcript_start"] + 400)

# Ensure BED order (start < end)
df["PAS_bed_start"] = df[["PAS_pos", "PAS_window"]].min(axis=1).astype(int)
df["PAS_bed_end"]   = df[["PAS_pos", "PAS_window"]].max(axis=1).astype(int)

# Set transcript_bed_start and transcript_bed_end strand-aware
df["transcript_bed_start"] = np.where(
    df["strand"] == "+",
    df[["transcript_end", "transcript_window"]].min(axis=1),
    df[["transcript_start", "transcript_window"]].min(axis=1)
).astype(int)

df["transcript_bed_end"] = np.where(
    df["strand"] == "+",
    df[["transcript_end", "transcript_window"]].max(axis=1),
    df[["transcript_start", "transcript_window"]].max(axis=1)
).astype(int)

# Prepare BED-format outputs
pas_df = df[["chrom", "PAS_bed_start", "PAS_bed_end", "strand"]]
endt_df = df[["chrom", "transcript_bed_start", "transcript_bed_end", "strand"]]

# Write to output files
pas_df.to_csv(pas_path, sep='\t', index=False, header=False)
endt_df.to_csv(endt_path, sep='\t', index=False, header=False)
