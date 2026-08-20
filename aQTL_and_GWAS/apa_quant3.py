import sys
import pandas as pd

counts_file = sys.argv[1]
bam_list = sys.argv[2]
raw_out = sys.argv[3]
usage_out = sys.argv[4]

# Load BAM names
with open(bam_list) as f:
    samples = [x.strip().split("/")[-1].replace(".bam","") for x in f]

# Read multicov output
df = pd.read_csv(counts_file, sep="\t", header=None)

chrom = df[0]
start = df[1]
end = df[2]
region = df[3]

counts = df.iloc[:,6:]
counts.columns = samples

# Build new RegionID
# preserve gene + coordinates
region_ids = (
    region + "|" +
    chrom + ":" +
    start.astype(str) + "-" +
    end.astype(str)
)

counts.index = region_ids

# Save raw count matrix
counts.to_csv(raw_out, sep="\t")

# Compute relative usage
genes = region.str.split("_").str[0]

usage = counts.groupby(genes.values).transform(
    lambda x: x / x.sum()
)

usage.index = region_ids

usage.to_csv(usage_out, sep="\t")
