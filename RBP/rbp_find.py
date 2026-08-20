#!/usr/bin/env python3

import os
import sys
from collections import defaultdict

# Usage:
# python map_rbp_fast.py ctrl_file eCLIPIDs eclip_dir out_file

ctrl_file = sys.argv[1]
eclip_ids_file = sys.argv[2]
eclip_dir = sys.argv[3]
out_file = sys.argv[4]

# chrom -> list of (begin, end, prot)
peaks_by_chrom = defaultdict(list)

print("Loading eCLIP files...", file=sys.stderr)

with open(eclip_ids_file) as f:
    clip_files = [line.strip() for line in f if line.strip()]

for name in clip_files:
    rbp_file = os.path.join(eclip_dir, name)
    prot = name.split("_")[0]

    with open(rbp_file) as f:
        for line in f:
            if not line.strip():
                continue

            fields = line.rstrip("\n").split()
            chrom = fields[0]
            begin = int(fields[1])
            end = int(fields[2])

            peaks_by_chrom[chrom].append((begin, end, prot))

print("Processing events...", file=sys.stderr)

with open(ctrl_file) as events, open(out_file, "w") as out:
    for event in events:
        event = event.rstrip("\n")
        if not event:
            continue

        fields = event.split("\t")

        rnaid = fields[0]
        brainid = fields[1]
        gene = fields[2]
        chrom = fields[3]
        pos = int(fields[4])
        nev = fields[5]
        signif = fields[6]
        tissue = fields[7]
        spec_stat = fields[8]
        statuss = fields[9]
        jobid = fields[10] if len(fields) > 10 else ""

        confidence = "Yes"

        for begin, end, prot in peaks_by_chrom.get(chrom, []):
            range_bot = begin - 1000
            range_top = end + 1000

            if begin <= pos <= end:
                distance = 0
                out.write(
                    f"{rnaid}\t{gene}\t{chrom}\t{pos}\t{tissue}\t{statuss}\t"
                    f"{confidence}\t{prot}\t{distance}\n"
                )

            elif range_bot <= pos <= range_top:
                distance = min(abs(begin - pos), abs(pos - end))
                out.write(
                    f"{rnaid}\t{gene}\t{chrom}\t{pos}\t{tissue}\t{statuss}\t"
                    f"{confidence}\t{prot}\t{distance}\n"
                )

print("Done.", file=sys.stderr)
