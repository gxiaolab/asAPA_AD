#!/usr/bin/python
import sys
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection
from statsmodels.stats.multitest import multipletests
##########

## python fdr_correction.py rbp_file_with_count_and_p_values
## Performs fdr correction on rbp p values and outputs adjusted p-values

##########

infile = sys.argv[1]

outfile = sys.argv[2]

df = pd.read_csv(infile, sep='\t')

p_values_ref = pd.to_numeric(df['Pvalue'])

corrected_p_values_ref = multipletests(p_values_ref, method='fdr_bh')[1]

df['adjusted_Pvalue'] = np.around(corrected_p_values_ref, 6)

df.to_csv(outfile, sep='\t', index=False)
