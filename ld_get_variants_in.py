import pandas as pd
import numpy as np
from collections import defaultdict
import sys
import re
from scipy.stats import ranksums
from time import time
import optparse
import os
import tabix

variant = sys.argv[1]

chrom, pos = variant.split()
pos = int(pos)

variants_in_ld = set()

primary_snp = f'{chrom}\t{pos}'
variants_in_ld.add(primary_snp)

for pop in ['EUR', 'AFR', 'EAS', 'SAS']:

    ld_file = f'./ld_database/file.txt.gz'  # Replace with the actual path to your LD file

    if not os.path.exists(ld_file):
        print("FILE DNE")
        continue

    tb = tabix.open(ld_file)

    try:
        for ld in tb.query(chrom, pos - 2, pos + 2):
            # chr12   11536   11537   1.0     1.0     .       11537:C:A       147241:T:C      +
            chrom, pos0, pos1, r2, dprime, dot, snp1, snp2, strand = ld

            if int(pos1) == pos and float(r2) >= 0.8 and float(dprime) >= 0.9:
                snp2_pos, snp2_ref, snp2_alt = snp2.split(':')
                snp2_pos = int(snp2_pos)
                snp2_id = f'{chrom}\t{snp2_pos}'
                variants_in_ld.add(snp2_id)
    except tabix.TabixError:
        continue

print('\n'.join(map(str, variants_in_ld)))
