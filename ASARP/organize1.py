#!/usr/bin/env python3

import re
import sys

if len(sys.argv) != 3:
    sys.stderr.write(f"Usage: {sys.argv[0]} <input_tmp2> <output_file>\n")
    sys.exit(1)

infile = sys.argv[1]
outfile = sys.argv[2]

with open(infile, "r") as f:
    text = f.read()

records = [r.strip() for r in text.split(";") if r.strip()]

front_pat = re.compile(
    r"""^TESTABLE:\s+
        (?P<signif>.+?)\s+candidate\ pair\ found!\s+
        (?P<gene>ENSG\d+(?:\.\d+)?):\s+
        (?P<gene_pos>\d+)\s+
        (?P<chrom>chr[\w]+)\s+
        (?P<gene_strand>[+-])\s+
        \((?P<event_block>.*)\)\s+
        VS\s+
        (?P<vs_pos>\d+):\s+
        (?P<rest>.*)$
    """,
    re.VERBOSE
)

at_pat = re.compile(r'\bAT:\s*3(?P<at_sign>[+-])\((?P<at_value>[^)]+)\)')
tail_pat = re.compile(
    r'^[^\s]+\s+\d+\s+(?P<snp>[ACGTN]+>[ACGTN]+)\s+nan\s+(?P<count1>\d+)\s+(?P<count2>\d+)$'
)

with open(outfile, "w") as out:
    for rec in records:
        rec_clean = " ".join(rec.split())

        m1 = front_pat.match(rec_clean)
        if not m1:
            sys.stderr.write(f"Could not parse front of record:\n{rec_clean}\n\n")
            continue

        m2 = at_pat.search(m1.group("event_block"))
        if not m2:
            sys.stderr.write(f"Could not parse AT field:\n{rec_clean}\n\n")
            continue

        m3 = tail_pat.match(m1.group("rest"))
        if not m3:
            sys.stderr.write(f"Could not parse tail of record:\n{rec_clean}\n\n")
            continue

        vals = [
            m1.group("signif"),
            m1.group("gene"),
            m1.group("gene_pos"),
            m1.group("chrom"),
            m2.group("at_sign"),
            m2.group("at_value"),
            m1.group("vs_pos"),
            m3.group("snp"),
            m3.group("count1"),
            m3.group("count2"),
        ]
        out.write("\t".join(vals) + "\n")
