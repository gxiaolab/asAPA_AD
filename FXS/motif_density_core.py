import argparse
import sys
import glob
from pysam import VariantFile
from pysam import FastaFile
import time
from time import strftime
from collections import defaultdict
import pandas as pd
import re
import math
import random
import numpy as np
from scipy.stats import kruskal, ks_2samp

# ======================
# Helper Functions
# ======================

fmr1 = ['ACUG', 'UGGA', 'ACUU', 'AGGA']

def rev_comp(seq):
    """Return the reverse complement of a DNA sequence."""
    comps = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(comps.get(base, base) for base in reversed(seq.upper()))

def dna_to_rna(seq):
    """Convert a DNA sequence to RNA by replacing T with U."""
    return seq.upper().replace('T', 'U')

def count_motifs(seq, motifs):
    """Count total motif occurrences in a sequence."""
    return sum(len(re.findall(motif, seq)) for motif in motifs)

def compute_gc_content(seq):
    gc = seq.count('G') + seq.count('C')
    return gc / len(seq) if len(seq) > 0 else 0

def downsample_until_balanced(seq_dict, alpha=0.05, bins=20, max_iter=1000):
    """
    Downsamples sequences in the largest groups to balance GC content across event categories.
    Only downsamples if GC content distributions differ significantly (p < alpha).

    Parameters:
        seq_dict (dict): {category: [sequences]}
        alpha (float): Significance level for Kruskal-Wallis test.
        bins (int): Not used directly but retained for consistency.
        max_iter (int): Max iterations to attempt balancing.

    Returns:
        dict: Balanced sequences per category.
    """
    # Compute initial GC content lists
    def get_gc_by_category(seq_dict):
        return {cat: [compute_gc_content(seq) for seq in seqs] for cat, seqs in seq_dict.items()}

    gc_content = get_gc_by_category(seq_dict)

    # Check if GC content differs significantly
    categories = list(gc_content.keys())
    values = list(gc_content.values())

    if len(categories) == 2:
        stat, p_value = ks_2samp(values[0], values[1])
    else:
        stat, p_value = kruskal(*values)

    if p_value >= alpha:
        return seq_dict  # Already balanced

    # Begin iterative downsampling
    seq_dict = {cat: seqs.copy() for cat, seqs in seq_dict.items()}

    for _ in range(max_iter):
        gc_content = get_gc_by_category(seq_dict)
        categories = list(gc_content.keys())
        values = list(gc_content.values())

        if len(categories) == 2:
            stat, p_value = ks_2samp(values[0], values[1])
        else:
            stat, p_value = kruskal(*values)

        if p_value >= alpha:
            break  # Balanced

        # Find category with highest number of sequences
        largest_cat = max(seq_dict.keys(), key=lambda k: len(seq_dict[k]))

        if len(seq_dict[largest_cat]) <= 1:
            break  # Avoid complete removal

        # Remove one random sequence from largest category
        seq_dict[largest_cat].pop(random.randint(0, len(seq_dict[largest_cat]) - 1))

    return seq_dict

def scan_motifs(sequence, event_category, window_size):
    """
    Scan sequence with sliding window and count motifs.
    Returns a list of (relative_position, count, category) tuples.
    """
    results = []
    seq_len = len(sequence)
    midpoint = (seq_len - window_size)// 2
    for i in range(seq_len - window_size + 1):
        window_seq = sequence[i:i + window_size]
        count = count_motifs(window_seq, fmr1)
        relative_pos = i - midpoint
        results.append((relative_pos, count, event_category))
    return results

def scan_motifs_percent_binned(sequence, event_category, motifs):
    results = []
    seq_len = len(sequence)
    bin_size = max(4, seq_len // 100)

    for bin_num in range(100):
        start = bin_num * bin_size
        end = (bin_num + 1) * bin_size if bin_num < 99 else seq_len
        bin_seq = sequence[start:end]
        motif_count = count_motifs(bin_seq, motifs)
        bin_len = end - start
        density = motif_count / bin_len if bin_len > 0 else 0
        results.append((bin_num + 1, density, event_category))
    return results

def scan_motifs_percent_binned_extended(sequence, event_category, motifs, core_len):
    """
    Similar to scan_motifs_percent_binned, but:
    - Supports percentage bins from -10% to 110%
    - Allows for later smoothing at the summary stage without clamping
    """
    results = []
    seq_len = len(sequence)
    bin_size = max(4, seq_len // 120)  # total 120 bins to support -10% to 110%

    for bin_num in range(-10, 111):
        bin_num_pseudo = bin_num + 10
        start = int((bin_num_pseudo / 100) * core_len)
        end = int(((bin_num_pseudo + 1) / 100) * core_len)

        # Make sure we're indexing valid regions
        start = max(0, start)
        end = min(seq_len, end)

        bin_seq = sequence[start:end]
        bin_len = len(bin_seq)
        motif_count = count_motifs(bin_seq, motifs)
        density = motif_count / bin_len if bin_len > 0 else 0

        results.append((bin_num, density, event_category))

    return results

def summarize_results_windowed_sliding(raw_data, window=10):
    """
    Summarize motif density using overlapping sliding windows.
    Each row represents the average and standard error of values within a window.

    Args:
        raw_data: list of (position, value, event_category)
        window: int, size of the sliding window

    Returns:
        summary: list of (midpoint_position, mean_density, event_category, standard_error)
    """
    # Organize by category
    grouped_by_category = defaultdict(list)
    for pos, val, cat in raw_data:
        grouped_by_category[cat].append((int(pos), val))

    summary = []

    for cat, data in grouped_by_category.items():
        data.sort(key=lambda x: x[0])
        positions = [p for p, _ in data]
        min_pos, max_pos = min(positions), max(positions)

        # Build a dictionary of values for fast lookup
        value_dict = defaultdict(list)
        for pos, val in data:
            value_dict[pos].append(val)

        for start in range(min_pos, max_pos - window + 2):  # inclusive of last window
            end = start + window
            window_vals = []
            for i in range(start, end):
                window_vals.extend(value_dict.get(i, []))
            n = len(window_vals)
            if n == 0:
                continue
            mean_val = sum(window_vals) / n
            se = (math.sqrt(sum((v - mean_val) ** 2 for v in window_vals) / (n - 1)) / math.sqrt(n)) if n > 1 else 0.0
            midpoint = start + window // 2
            summary.append((midpoint, mean_val, cat, se))

    summary.sort(key=lambda x: (x[2], x[0]))  # sort by category then position
    return summary

def write_summary(filename, summary_data, position_label):
    with open(filename, 'w') as f:
        f.write(f"{position_label}\tMean_Motif_Density\tEvent_Category\tStandard_Error\n")
        for row in summary_data:
            f.write("\t".join(str(val) for val in row) + "\n")

# ======================
# Main Script
# ======================

def main():
    parser = argparse.ArgumentParser(description='Extract RNA sequences and scan for motifs around PAS and transcript sites.')
    parser.add_argument('-i', '--input', required=True, help='Input file with PAS and transcript coordinates')
    parser.add_argument('-g', '--genome', required=True, help='Path to genome FASTA file')
    parser.add_argument('-o', '--output_prefix', required=True, help='Output file prefix')
    parser.add_argument('-w', '--window', type=int, default=300, help='Window size around PAS and transcript sites')
    parser.add_argument('--sliding', type=int, required=True, help='Sliding window size (must be between 4 and 600)')
    parser.add_argument('--smooth', type=int, default=10, help='Summary Smoothing window size')
    args = parser.parse_args()

    if not (4 <= args.sliding <= 600):
        sys.exit("Sliding window size must be between 4 and 600.")

    # Load input data
    df = pd.read_csv(args.input, sep='\t')

    genome = FastaFile(args.genome)

    all_pas_results = []
    all_tx_results = []
    all_span_results = []

    all_pas_seqs = {'Shortened': [], 'Lengthened': [], 'Not': []}
    all_tx_seqs = {'Shortened': [], 'Lengthened': [], 'Not': []}
    all_span_seqs = {'Shortened': [], 'Lengthened': [], 'Not': []}

    processed_span = set()

    for idx, row in df.iterrows():
        chrom = row['chrom']
        strand = row['strand']
        event_category = row['event_category']
        
        #pas_pos = int(row['region_start'])
        #tx_pos = int(row['region_end'])

        pas_pos = int(row['PAS'])
        tx_pos = int(row['txend'])

        #pas_pos = int(row['PAS_start']) if strand == '-' else int(row['PAS_end'])
        #tx_pos = int(row['transcript_start']) if strand == '-' else int(row['transcript_end'])

        span_key = (chrom, strand, event_category, pas_pos, tx_pos)

        if span_key in processed_span:
            continue

        skip_span = span_key in processed_span

        if not skip_span:
            processed_span.add(span_key)

        span_start = min(pas_pos, tx_pos)
        span_end = max(pas_pos, tx_pos)
        #span_start = pas_pos
        #span_end = tx_pos

        span_length = span_end - span_start
     
        pad = int(0.1 * span_length)

        expanded_start = span_start - pad
        expanded_end = span_end + pad

        if strand == '-':
            span_seq = rev_comp(genome.fetch(chrom, expanded_start, expanded_end))
            span_rna = dna_to_rna(span_seq)
            all_span_seqs[event_category].append(span_rna)
        else:
            span_seq = genome.fetch(chrom, expanded_start, expanded_end)
            span_rna = dna_to_rna(span_seq)
            all_span_seqs[event_category].append(span_rna)

    # Apply GC content correction only if distributions differ significantly
    balanced_span = downsample_until_balanced(all_span_seqs, alpha=0.05, bins=20)

    for cat, seqs in balanced_span.items():
        for seq in seqs:
            span_len = len(seq) * (100/120)  # Estimate original span_len
            all_span_results.extend(scan_motifs_percent_binned_extended(seq, cat, fmr1, span_len))

    # Summarize and output
    span_summary = summarize_results_windowed_sliding(all_span_results, 10)

    write_summary(f"{args.output_prefix}_SPAN_summary.tsv", span_summary, "Percent_Bin")

if __name__ == "__main__":
    main()


