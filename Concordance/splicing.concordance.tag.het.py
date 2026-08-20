#!/usr/bin/python

import sys
import argparse
import glob
from time import strftime
import os
import time

#import GenomeFetch as gf
#gf = gf.GenomeFetch('hg19')
import subprocess

###########
#- inputs:
#	- list of tag snps (bed format)
#		*info column: refcount:varcount[:othercount]|chrm:exon_st1:exon_end:strd|genotypeSource
#	- list of gt snps (bed format)
#		*info column: gt => 0/0; 0/1; 1/1
#- output:
#	- Si for each input snps
#
# To reduce some pairwise comparisons b/t SNPs:
#	- Tag SNP + candidate SNPs w/in 500bp near exon-intron boundaries
#	- if "INF" for -s option, use anno => only search for snp pairs w/in the same gene
#
#	there are complicated GT such like chr1.16388709. In this case, filter out in splicing.concordance.py!
#----------
# NEW: use genotype source info from tag snv.
#	=> if source is from GT, max d is as expected @ 0.5
#	=> if source is from db, max d depends on AR thresholds used to determine trustable heterozygous ASAS tag snvs.
###########

parser = argparse.ArgumentParser(description='Script descriptions here')
parser.add_argument('-i', metavar='inf', required=True, help='Input causal candidate bed file')
parser.add_argument('-d', metavar='indiv', required=True, help='Input individual ID')
parser.add_argument('-t', metavar='tag', required=True, help='tag snp bed file')
parser.add_argument('-m', metavar='maxD', required=True, help='max d for RNA-seq defined tag snvs')
parser.add_argument('-o', metavar='outf', required=True, help='Output file')
parser.add_argument('-a', metavar='anno', required=True, help='gene anno bed file')
parser.add_argument('-s', metavar='search', required=True, help='max dist in nt from candidate causal snp to the AS exon to be tested; input "INF" to test all possible snp pairs w/in the same gene')
parser.add_argument('-b', metavar='hettag', required=True, help='tag snps that are heterozygous in the individual')

opts = parser.parse_args()
print('Inf candidate causal snp: %s' % opts.i)
print('indiv id: %s' % opts.d)
print('tag snp: %s' % opts.t)
print('max d for RNA-seq defined tag snvs: %s' % opts.m)
print('anno: %s' % opts.a)
print('search: %s' % opts.s)
print('Outf: %s' % opts.o)
print('het tag snps: %s' % opts.b)

maxD = float(opts.m)

try:
	maxdist = int(opts.s)
except ValueError:
	maxdist = 'INF'

def snpInGene(inf, anno):
    res = {}

    p = subprocess.run(
        ['intersectBed',
         '-wo', '-a', inf, '-b', anno],
        capture_output=True
    )

    stdout = p.stdout.decode()
    stderr = p.stderr.decode()

    # Keep these critical print statements
    #print(stdout, stderr)

    if stderr:
        print(1, inf)
        sys.exit()

    # New guard: no overlaps
    if not stdout.strip():
        print(f"WARNING: no overlaps found between {inf} and {anno}")
        return res

    for line in stdout.strip().split('\n'):
        l = line.split('\t')

        if len(l) < 11:
            print(f"WARNING: malformed intersect line skipped: {line}")
            continue

        g = l[-4].split('|')[0]

        if g not in res:
            res[g] = {}

        res[g][tuple(l[:-7] + [l[-2]])] = 1

    return res

def readHetTags(hettag):
    res = {}
    with open(hettag, 'r') as f:
        for line in f:
            chr, start, end, gt = line.strip().split('\t')
            if chr not in res:
                res[chr] = [end]
            else:
                res[chr].append(end)
    return res

def main(argv):
	print('job starts', strftime('%a, %d %b %Y %I:%M:%S'))
	start_time = time.time()
 
 	# list of heterozygous tag snps
	hettag = readHetTags(opts.b)

	#always search snp pairs w/in genes
	#candidate causal snp overlap w/ gene anno
	candidate = snpInGene(opts.i,opts.a) #candidate[gene]:snp-bed-file
	#tag snp overlap w/ gene anno
	tag = snpInGene(opts.t,opts.a)
 


	out = open(opts.o,'w')
	out.write('causalCandidate\tgenotype\texon\tindiv|tagGTsource\ttagSNP\tcoverage\tref\tvar\tallelicR\tdi\tdi2\tSi\n')
	if maxdist == 'INF': 
		for g in set(candidate.keys()).intersection(set(tag.keys())): #intersection means candidate and tag are in the same gene
			for cl in candidate[g]:
				csnp = '{}.{}.{}'.format(cl[0],cl[2],cl[-1]) #with strd ==> check for same strd later in get.causal.py!
				gt = cl[3]
				try:
					gtreal = gt.split('|')[0]
				except IndexError: print('index err',gt); continue
				for tl in tag[g]:
					tsnp = '{}.{}.{}'.format(tl[0],tl[2],tl[-1]) #with strd ==> check for same strd later in get.causal.py!
					counts,exon,source = tl[3].split('|')
					chrm,st1,end,strd = exon.split(':')
					if cl[-1] == tl[-1] and cl[-1] == strd:
						counts = list(map(int,counts.split(':')))
						ref,var = counts[:2]
						tot = sum(counts)
						ri = 1.*ref/tot
						di = abs(0.5-ri)
						if tsnp == csnp:# or gtreal in ('0/1', '1/0'):
							if source == 'GT': si = di**2/0.5**2
							else: si = di**2/(maxD)**2
						#elif tl[2] in hettag[chrm]:
						elif tl[2] in hettag.get(chrm, []):
							if gtreal in ('0/1', '1/0'):
								if source == 'GT': si = di**2/0.5**2
								else: si = di**2/(maxD)**2
							elif gtreal in ('0/0', '1/1'):
								if source == 'GT': si = 1-di**2/0.5**2
								else: si = 1-di**2/(maxD)**2
						elif gtreal in ('0/0', '1/1', '0/1', '1/0'): continue
						else: print('other gt',csnp,gt); continue
						out.write('{}\t{}\t{}\t{}|{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(csnp,gt,exon,opts.d,source,tsnp,tot,ref,var,ri,di,di**2,si))
	else:
		for g in set(candidate.keys()).intersection(set(tag.keys())): #intersection means candidate and tag are in the same gene
			for cl in candidate[g]:
				csnp = '{}.{}.{}'.format(cl[0],cl[2],cl[-1]) #with strd ==> check for same strd later in get.causal.py!
				gt = cl[3]
				try:
					gtreal = gt.split('|')[0]
				except IndexError: print('index err',gt); continue
				for tl in tag[g]:
					tsnp = '{}.{}.{}'.format(tl[0],tl[2],tl[-1]) #with strd ==> check for same strd later in get.causal.py!  
					counts,exon,source = tl[3].split('|')
					chrm,st1,end,strd = exon.split(':')
					if cl[-1] == tl[-1] and cl[-1] == strd:
						st1,end,pos = int(st1), int(end), int(cl[2])
						if abs(st1-pos) <= maxdist or abs(end-pos) <= maxdist or pos >= st1 and pos <= end:
							counts = map(int,counts.split(':'))
							ref,var = counts[:2]
							tot = sum(counts)
							ri = 1.*ref/tot
							di = abs(0.5-ri)
							if tsnp == csnp or gtreal in ('0/1', '1/0'): 
								if source == 'GT': si = di**2/0.5**2
								else: si = di**2/(maxD)**2
							elif gtreal in ('0/0', '1/1'):
								if source == 'GT': si = 1-di**2/0.5**2
								else: si = 1-di**2/(maxD)**2
							#1/2/2018: haven't tried the following 2 lines yet...
							#if tsnp == csnp or gtreal[0] != gtreal[-1]: si = di**2/0.5**2
							#elif gtreal[0] == gtreal[-1]: si = 1-di**2/0.5**2
							else: print('other gt',csnp,gt); continue
							out.write('{}\t{}\t{}\t{}|{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(csnp,gt,exon,opts.d,source,tsnp,tot,ref,var,ri,di,di**2,si))

	print("--- %s seconds ---" % (time.time() - start_time))
	print('DONE!', strftime('%a, %d %b %Y %I:%M:%S'))

if __name__ == '__main__':
	main(sys.argv[1:])
