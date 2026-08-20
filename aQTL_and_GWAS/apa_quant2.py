import sys
from collections import defaultdict

bed_file = sys.argv[1]
annotation = sys.argv[2]
out_bed = sys.argv[3]

# -----------------------------
# cluster function
# -----------------------------
def cluster_sites(sites, window=50):

    sites = sorted(set(sites))
    clusters = []
    current = [sites[0]]

    for s in sites[1:]:
        if s - current[-1] <= window:
            current.append(s)
        else:
            clusters.append(current)
            current = [s]

    clusters.append(current)

    return [int(sum(c)/len(c)) for c in clusters]


# -----------------------------
# build annotation lookup
# -----------------------------
plus_lookup = {}
minus_lookup = {}

na_counter = 1

with open(annotation, encoding="utf-8", errors="replace") as f:

    for line in f:

        if line.startswith("original"):
            continue

        L = line.strip().split()

        if len(L) < 9:
            continue

        chrom = L[0]
        pos = int(L[1])
        strand = L[2]

        # ----- FORMAT WITH GENE -----
        if L[3].startswith("ENSG"):

            gene = L[3]
            rna_start = int(float(L[8]))
            rna_end = int(float(L[9]))

        # ----- FORMAT WITHOUT GENE -----
        elif L[3] == "na":
            
            gene = f"NA{na_counter}"
            na_counter += 1

            rna_start = int(float(L[8]))
            rna_end = int(float(L[9]))

        elif L[3] == "NO":
            
            gene = L[5]
            na_counter += 1

            rna_start = int(float(L[8]))
            rna_end = int(float(L[9]))

        else:

            gene = f"NA{na_counter}"
            na_counter += 1

            rna_start = int(float(L[5]))
            rna_end = int(float(L[6]))

        if strand == "+":
            key = (chrom, pos, rna_end, strand)
            plus_lookup[key] = gene

        else:
            key = (chrom, pos, rna_start, strand)
            minus_lookup[key] = gene

# collect PAS per gene
gene_sites = defaultdict(list)
gene_chrom = {}
gene_strand = {}

with open(bed_file) as f:

    for line in f:

        chrom, c2, c3, c4, c5, strand = line.strip().split()

        c2 = int(c2)
        c3 = int(c3)
        c4 = int(c4)
        c5 = int(c5)

        if strand == "+":
            key = (chrom, c3, c5, strand)
            gene = plus_lookup.get(key)

            site = c5

        else:
            key = (chrom, c2, c4, strand)
            gene = minus_lookup.get(key)

            site = c2

        if gene is None:
            continue

        gene_sites[gene].append(site)
        gene_chrom[gene] = chrom
        gene_strand[gene] = strand

# build UTR segments
out = open(out_bed,"w")

for gene, sites in gene_sites.items():

    clusters = cluster_sites(sites)

    if len(clusters) < 2:
        continue

    chrom = gene_chrom[gene]
    strand = gene_strand[gene]

    if strand == "+":
        ordered = sorted(clusters)
    else:
        ordered = sorted(clusters, reverse=True)

    segments = []

    for i in range(len(ordered)-1):

        a = ordered[i]
        b = ordered[i+1]

        if strand == "+":
            start = a
            end = b
        else:
            start = b
            end = a

        if abs(end-start) < 50:
            continue

        region_id = f"{gene}_UTR{i+1}"

        out.write(
            "\t".join(map(str,[
                chrom,
                start,
                end,
                region_id,
                0,
                strand
            ]))+"\n"
        )

out.close()
