import sys
from collections import defaultdict
from optparse import OptionParser


def cluster_sites(sites, window=50):

    if len(sites) == 0:
        return []

    sites = sorted(sites)
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


if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option("-i", dest="annotation")
    parser.add_option("-b", dest="bed")
    parser.add_option("-o", dest="events")

    (options, args) = parser.parse_args()

    if not options.annotation or not options.bed or not options.events:
        parser.error("Usage: -i annotation -b bed_output -o event_output")

    annotation = options.annotation

    events_file = open(options.events, "w")
    bed_file = open(options.bed, "w")

    gene_transcripts = defaultdict(list)
    transcripts = {}

    print("Parsing annotation...")

    with open(annotation) as f:

        for line in f:

            if line.startswith("#"):
                continue

            L = line.strip().split()

            if len(L) < 11:
                continue

            transcript_id, chrom, strand, tx_start, tx_end, cds_start, cds_end, exon_n, exon_starts, exon_ends, gene_id = L[:11]

            exon_starts = [int(x) for x in exon_starts.split(",") if x]
            exon_ends   = [int(x) for x in exon_ends.split(",") if x]

            if len(exon_starts) == 0:
                continue

            transcripts[transcript_id] = {
                "gene": gene_id,
                "chrom": chrom,
                "strand": strand,
                "starts": exon_starts,
                "ends": exon_ends
            }

            gene_transcripts[gene_id].append(transcript_id)

    print("Genes loaded:", len(gene_transcripts))

    gene_counter = 0
    region_counter = 0

    for gene_id, tx_list in gene_transcripts.items():

        chrom = transcripts[tx_list[0]]["chrom"]
        strand = transcripts[tx_list[0]]["strand"]

        polyA_sites = []
        upstream_candidates = []

        terminal_exons = []
        start_support = defaultdict(int)

        for tx in tx_list:

            T = transcripts[tx]

            starts = T["starts"]
            ends = T["ends"]

            if len(starts) < 1:
                continue

            if strand == "+":
                start = starts[-1]
                end = ends[-1]
                prev_end = ends[-2] if len(ends) > 1 else start

            else:
                start = starts[0]
                end = ends[0]
                prev_end = starts[1] if len(starts) > 1 else end

            exon_length = abs(end - start)

            # FILTER 1: remove very short terminal exons
            if exon_length < 50:
                continue

            terminal_exons.append((start, end, prev_end))

            start_support[start] += 1

            polyA_sites.append(end if strand == "+" else start)

            upstream_candidates.append(prev_end)

        if len(terminal_exons) == 0:
            continue

        # FILTER 2: remove singleton terminal exon starts
        filtered_exons = []

        for start, end, prev_end in terminal_exons:

            if start_support[start] < 2:
                continue

            filtered_exons.append((start, end, prev_end))

        if len(filtered_exons) == 0:
            continue

        polyA_sites = list(set(polyA_sites))
        clustered_sites = cluster_sites(polyA_sites)

        if len(clustered_sites) < 2:
            continue

        if strand == "+":
            ordered = sorted(clustered_sites)
            upstream = min(upstream_candidates)
        else:
            ordered = sorted(clustered_sites, reverse=True)
            upstream = max(upstream_candidates)

        gene_counter += 1

        segments = []

        if strand == "+":
            segments.append((upstream, ordered[0]))
        else:
            segments.append((ordered[0], upstream))

        for i in range(len(ordered)-1):

            a = ordered[i]
            b = ordered[i+1]

            if strand == "+":
                segments.append((a,b))
            else:
                segments.append((b,a))

        for i,(start,end) in enumerate(segments):

            if abs(end-start) < 50:
                continue

            region_id = f"{gene_id}_UTR{i+1}"

            events_file.write(
                "\t".join(map(str,[
                    "APA_UTR",
                    chrom,
                    gene_id,
                    region_id,
                    start,
                    end,
                    strand
                ]))+"\n"
            )

            bed_file.write(
                "\t".join(map(str,[
                    chrom,
                    start,
                    end,
                    region_id,
                    0,
                    strand
                ]))+"\n"
            )

            region_counter += 1

        # ---------- ALE detection ----------

        start_positions = [x[0] for x in filtered_exons]
        start_clusters = cluster_sites(list(set(start_positions)), window=50)

        if len(start_clusters) > 1:

            if strand == "+":
                ordered_starts = sorted(start_clusters)
            else:
                ordered_starts = sorted(start_clusters, reverse=True)

            for i in range(len(ordered_starts)-1):

                a = ordered_starts[i]
                b = ordered_starts[i+1]

                if strand == "+":
                    start = a
                    end = b
                else:
                    start = b
                    end = a

                if abs(end-start) < 50:
                    continue

                region_id = f"{gene_id}_ALE{i+1}"

                events_file.write(
                    "\t".join(map(str,[
                        "APA_ALE",
                        chrom,
                        gene_id,
                        region_id,
                        start,
                        end,
                        strand
                    ]))+"\n"
                )

                bed_file.write(
                    "\t".join(map(str,[
                        chrom,
                        start,
                        end,
                        region_id,
                        0,
                        strand
                    ]))+"\n"
                )

                region_counter += 1

    events_file.close()
    bed_file.close()

    print("APA genes detected:", gene_counter)
    print("UTR/ALE segments written:", region_counter)
