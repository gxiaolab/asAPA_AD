import pandas as pd
import numpy as np
import math
import sys

fxtas = sys.argv[1]
neurobiobank = sys.argv[2]
out_file = sys.argv[3]

# Load datasets
fx = pd.read_csv(fxtas, sep='\t')
nb = pd.read_csv(neurobiobank, sep='\t')

# Create a unique event ID
fx['event_id'] = fx['chr'].astype(str) + "_" + fx['PAS_end'].astype(str) + "_" + fx['transcript_end'].astype(str)
nb['event_id'] = nb['chr'].astype(str) + "_" + nb['PAS_end'].astype(str) + "_" + nb['transcript_end'].astype(str)

# Subset necessary columns
fx_sub = fx[['event_id', 'chr', 'PAS_start', 'PAS_end', 'transcript_start', 'transcript_end', 'strand',
             'LFC_Difference', 'Event_Category']]
nb_sub = nb[['event_id', 'chr', 'PAS_start', 'PAS_end', 'transcript_start', 'transcript_end', 'strand',
             'LFC_Difference', 'Event_Category']]

# Merge datasets
merged = pd.merge(fx_sub, nb_sub, on='event_id', how='outer', suffixes=('_FX', '_NB'))

# Function to determine final event category
def consolidate_category(row):
    fx_cat = row['Event_Category_FX']
    nb_cat = row['Event_Category_NB']
    
    if fx_cat in ['Shortened', 'Lengthened'] and nb_cat in ['Shortened', 'Lengthened']:
        if fx_cat == nb_cat:
            print(f"{fx_cat} - Agree")
            return fx_cat
        else:
            return 'Disagree'
    elif fx_cat in ['Shortened', 'Lengthened']:
        return fx_cat
    elif nb_cat in ['Shortened', 'Lengthened']:
        return nb_cat
    elif fx_cat == 'Not' and nb_cat == 'Not':
        return 'Not'
    elif (fx_cat == 'Not' and pd.isna(nb_cat)) or (pd.isna(fx_cat) and nb_cat == 'Not'):
        return 'Not'
    elif pd.isna(fx_cat) and pd.isna(nb_cat):
        return 'NA'
    else:
        return 'ALL_FAIL'

# Function to determine final Delta_usage
def consolidate_delta(row):
    fx_delta = row['LFC_Difference_FX']
    nb_delta = row['LFC_Difference_NB']
    fx_cat = str(row['Event_Category_FX'])
    nb_cat = str(row['Event_Category_NB'])
    final_cat = row['Final_Category']

    # If both files agree (Shortened/Shortened, Lengthened/Lengthened, Not/Not)
    if fx_cat == nb_cat and final_cat in ['Shortened', 'Lengthened', 'Not']:
        if pd.notna(fx_delta) and pd.notna(nb_delta):
            return (fx_delta + nb_delta) / 2
        elif pd.notna(fx_delta):
            return fx_delta
        elif pd.notna(nb_delta):
            return nb_delta

    # If one is Shortened/Lengthened and the other is Not or NA → take the Shortened/Lengthened value
    if final_cat in ['Shortened', 'Lengthened']:
        if fx_cat == final_cat and pd.notna(fx_delta):
            return fx_delta
        elif nb_cat == final_cat and pd.notna(nb_delta):
            return nb_delta

    # If one is Not and the other is NA → take the Not value
    if final_cat == 'Not':
        if fx_cat == 'Not' and pd.notna(fx_delta):
            return fx_delta
        elif nb_cat == 'Not' and pd.notna(nb_delta):
            return nb_delta

    # If they Disagree → return NA
    return 'NA'


# Apply logic
merged['Final_Category'] = merged.apply(consolidate_category, axis=1)
merged['Final_Delta'] = merged.apply(consolidate_delta, axis=1)

# Prepare final output
output = pd.DataFrame({
    'chr': merged['chr_FX'],
    'PAS_start': merged['PAS_start_FX'],
    'PAS_end': merged['PAS_end_FX'],
    'transcript_start': merged['transcript_start_FX'],
    'transcript_end': merged['transcript_end_FX'],
    'strand': merged['strand_FX'],
    'PAS_coverage_FXS': 0,
    'PAS_coverage_Control': 0,
    'Transcript_coverage_FXS': 0,
    'Transcript_coverage_Control': 0,
    'LFC_PAS': 0,
    'LFC_Transcript': 0,
    'LFC_Difference': merged['Final_Delta'],
    'Event_Category': merged['Final_Category']
})

# Save to output file
output.to_csv(out_file, sep='\t', index=False)

print("Finished. Output saved")

