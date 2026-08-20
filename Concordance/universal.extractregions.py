#!/usr/bin/python

import sys

file = sys.argv[1]
out = sys.argv[2]

with open(file, 'r') as f:
    with open(out, 'w') as outfile:
        for line in f:
            chr, pos, strd, gene, region, type, ref, alt = line.strip().split('\t')
            start, end = region.split('-')
            outfile.write(f"{chr}\t{pos}\t{strd}\t{gene}\t{start}:{end}\t{type}\t{ref}\t{alt}\n")
        
    
