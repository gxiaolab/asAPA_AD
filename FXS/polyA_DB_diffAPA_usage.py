#!/usr/bin/python
import pandas as pd
import numpy as np
import sys

# Input arguments
reference_file = sys.argv[1]  # reference.bed
FXS_PAS_file = sys.argv[2]  # FXS_avg_PAS_coverage.bed
CTRL_PAS_file = sys.argv[3]  # Control_avg_PAS_coverage.bed
FXS_transcript_file = sys.argv[4]  # FXS_avg_transcript_coverage.bed
CTRL_transcript_file = sys.argv[5]  # Control_avg_transcript_coverage.bed
output_file = sys.argv[6]  # Final output file

# Load reference file
reference = pd.read_csv(reference_file, sep='\t', header=None, names=["chr", "PAS_start", "PAS_end", "transcript_start", "transcript_end", "strand"])
print(reference["transcript_start"].dtype)
print(reference["transcript_end"].dtype)

# Load coverage data
FXS_PAS = pd.read_csv(FXS_PAS_file, sep='\t', header=None, names=["chr", "start", "end", "strand", "PAS_coverage_FXS"])
CTRL_PAS = pd.read_csv(CTRL_PAS_file, sep='\t', header=None, names=["chr", "start", "end", "strand", "PAS_coverage_Control"])
FXS_transcript = pd.read_csv(FXS_transcript_file, sep='\t', header=None, names=["chr", "start", "end", "strand", "Transcript_coverage_FXS"])
CTRL_transcript = pd.read_csv(CTRL_transcript_file, sep='\t', header=None, names=["chr", "start", "end", "strand", "Transcript_coverage_Control"])

# Convert "." to 0 and ensure numeric coverage
for df in [FXS_PAS, CTRL_PAS, FXS_transcript, CTRL_transcript]:
    df.iloc[:, -1] = pd.to_numeric(df.iloc[:, -1], errors='coerce').fillna(0)

# Function to match coverage based on PAS and transcript coordinates
def match_coverage(reference, PAS_FXS, PAS_CTRL, transcript_FXS, transcript_CTRL):
    print(transcript_FXS["start"].dtype)
    print(transcript_FXS["end"].dtype)
    print(transcript_CTRL["start"].dtype)
    print(transcript_CTRL["end"].dtype)
    merged = reference \
        .merge(PAS_FXS, left_on=["chr", "PAS_start", "PAS_end", "strand"], right_on=["chr", "start", "end", "strand"], how="left") \
        .drop(columns=["start", "end"]) \
        .merge(PAS_CTRL, left_on=["chr", "PAS_start", "PAS_end", "strand"], right_on=["chr", "start", "end", "strand"], how="left") \
        .drop(columns=["start", "end"]) \
        .merge(transcript_FXS, left_on=["chr", "transcript_start", "transcript_end", "strand"], right_on=["chr", "start", "end", "strand"], how="left") \
        .drop(columns=["start", "end"]) \
        .merge(transcript_CTRL, left_on=["chr", "transcript_start", "transcript_end", "strand"], right_on=["chr", "start", "end", "strand"], how="left") \
        .drop(columns=["start", "end"])

    print(merged.head())
    return merged.fillna(0)  # Replace NaN with 0

# Merge all coverage data
merged_df = match_coverage(reference, FXS_PAS, CTRL_PAS, FXS_transcript, CTRL_transcript)

# Function to compute log2 fold change (LFC)
def compute_LFC(FXS, Control):
    """ Compute log2 fold change (LFC = log2(FXS / Control)), requiring a minimum coverage of 1. """
    if FXS < 5 or Control < 5:
        return "NA"  # If coverage is too low, return NA
    return np.log2(FXS / Control)

# Compute LFC for PAS and transcript sites
merged_df["LFC_PAS"] = merged_df.apply(lambda row: compute_LFC(row["PAS_coverage_FXS"], row["PAS_coverage_Control"]), axis=1)
merged_df["LFC_Transcript"] = merged_df.apply(lambda row: compute_LFC(row["Transcript_coverage_FXS"], row["Transcript_coverage_Control"]), axis=1)

# Compute difference between PAS and transcript LFC
def compute_LFC_difference(row):
    if row["LFC_PAS"] == "NA" or row["LFC_Transcript"] == "NA":
        return "NA"
    return row["LFC_PAS"] - row["LFC_Transcript"]

merged_df["LFC_Difference"] = merged_df.apply(compute_LFC_difference, axis=1)

# Categorize transcript event
# CHANGE CUTOFFS HERE
# First run was 0.5, second run was 1, third was 0.585
def categorize_event(LFC_diff):
    if LFC_diff == "NA":
        return "NA"
    # LFC_diff > (np.log2(1.5))
    elif LFC_diff > 0.5:
        return "Shortened"
    elif LFC_diff < -0.5:
        return "Lengthened"
    else:
        return "Not"

merged_df["Event_Category"] = merged_df["LFC_Difference"].apply(categorize_event)

# Save the final output
merged_df.to_csv(output_file, sep='\t', index=False)

