#!/usr/bin/python
import pandas as pd
import glob
import sys

file_path = sys.argv[1]
out_path = sys.argv[2]

# Function to compute average coverage
def average_coverage(file_pattern):
    files = glob.glob(file_pattern)  # Get all matching files
    df_list = []
    
    for f in files:
        df = pd.read_csv(f, sep='\t', header=None, names=["chr", "start", "end", "strand", "coverage"])
        
        # Convert "." to 0 in coverage column
        df["coverage"] = pd.to_numeric(df["coverage"], errors='coerce').fillna(0)
        
        df_list.append(df)

    # Merge on genomic coordinates and take mean coverage across files
    merged_df = pd.concat(df_list).groupby(["chr", "start", "end", "strand"], as_index=False)["coverage"].mean()
    return merged_df

# Compute averages for PAS and transcript sites for both conditions
FX_PAS_avg = average_coverage(f"{file_path}/*_FXS_PAS_coverage.bed")
CTRL_PAS_avg = average_coverage(f"{file_path}/*_Control_PAS_coverage.bed")
FX_transcript_avg = average_coverage(f"{file_path}/*_FXS_transcript_coverage.bed")
CTRL_transcript_avg = average_coverage(f"{file_path}/*_Control_transcript_coverage.bed")

# Save results
FX_PAS_avg.to_csv(f"{out_path}/FXS_avg_PAS_coverage.bed", sep='\t', index=False, header=False)
CTRL_PAS_avg.to_csv(f"{out_path}/Control_avg_PAS_coverage.bed", sep='\t', index=False, header=False)
FX_transcript_avg.to_csv(f"{out_path}/FXS_avg_transcript_coverage.bed", sep='\t', index=False, header=False)
CTRL_transcript_avg.to_csv(f"{out_path}/Control_avg_transcript_coverage.bed", sep='\t', index=False, header=False)

