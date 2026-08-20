import sys

universal_tag = sys.argv[1]
rdd_file = sys.argv[2]
outfile = sys.argv[3]

# read file from RDD 
snps = {}
with open(rdd_file, 'r') as f:
    for line in f:
        chr, pos, alleles, x, counts, y = line.split(' ')
        ref_count, var_count, var2_count = counts.split(':')
        snps[(chr, str(int(pos)-1), alleles)] = [ref_count, var_count, var2_count]
        
with open(outfile, 'w') as out:
    with open(universal_tag, 'r') as f:
        for line in f:
            if line.startswith('chrm'):
                continue
            chr, pos, strnd, gene, region, type, ref, alt = line.strip().split('\t')
            if (chr, pos, '{}>{}'.format(ref, alt)) not in snps: continue
            ref_count, var_count, var2_count = snps[(chr, pos, '{}>{}'.format(ref, alt))]
            col4 = '{}:{}:{}|{}:{}:{}:{}|GT'.format(str(ref_count),str(var_count),str(var2_count),chr, region.split(':')[0], region.split(':')[1], strnd)
            out.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(chr, pos, str(int(pos)+1), col4, str(0), strnd))
