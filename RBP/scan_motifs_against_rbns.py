#!/usr/bin/python

import argparse
import sys
import glob
import time
import itertools
from time import strftime
from collections import defaultdict as dd

def get_bases_degenerate(base):
	nt_codes = {
		'W': ['A', 'T'],
		'S': ['C', 'G'],
		'M': ['A', 'C'],
		'K': ['G', 'T'],
		'R': ['A', 'G'],
		'Y': ['C', 'T'],
		'B': ['C', 'G', 'T'],
		'D': ['A', 'G', 'T'],
		'H': ['A', 'C', 'T'],
		'V': ['A', 'C', 'G'],
		'N': ['A', 'C', 'G', 'T']
	}
	return nt_codes[base]

def main(argv):
	print('job starts', strftime('%a, %d %b %Y %I:%M:%S'))
	start_time = time.time()

	#read RBNS motifs into a dictionary
	rbp_motifs = dd(list)
	with open(opts.r) as f:
		rbps = []
		for line in f:
			if line.startswith('A1CF'):
				line = line.rstrip().split('\t')
				for i in range(0, len(line), 3):
					rbps.append(line[i].split()[0])
			else:
				line = line.rstrip().split('\t')
				for i in range(0, len(line), 3):
					rbp = rbps[i//3]
					motif = line[i]
					if not motif: continue
					rbp_motifs[rbp].append(motif)

	# print(rbp_motifs)

	with open(opts.o, 'w') as out:
		out.write('Motif\tDegenerate_Motif\tRBP\n')
		for motif_file in glob.glob('{}/*'.format(opts.i)):
			motif = motif_file.split('/')[-1].split('.')[0]
			dna_bases = 'ACGT'
			if all(base in dna_bases for base in motif): #Check whether motif has degenrate bases
				for rbp in rbp_motifs:
					if motif in rbp_motifs[rbp]:
						# print(rbp, motif)
						out.write('{}\t{}\t{}\n'.format(motif, 'NA', rbp))
			else:
				deg_idx = [i for i in range(len(motif)) if not motif[i] in dna_bases]
				if len(deg_idx) == 1:
					deg_bases = get_bases_degenerate(motif[deg_idx[0]])
					for b in deg_bases:
						deg_motif = motif.replace(motif[deg_idx[0]], b)
						#print(deg_motif)
						for rbp in rbp_motifs:
							if deg_motif in rbp_motifs[rbp]:
								#print(rbp, motif, deg_motif)
								out.write('{}\t{}\t{}\n'.format(motif, deg_motif, rbp))
				else:
					iter_items = []
					for idx in deg_idx:
						iter_items.append(get_bases_degenerate(motif[idx]))
					deg_base_combos = list(itertools.product(*iter_items))
					for c in deg_base_combos:
						deg_motif = motif
						for i in range(len(deg_idx)):
							deg_motif = deg_motif[:deg_idx[i]] + c[i] + deg_motif[deg_idx[i]+1:]
						if deg_motif in rbp_motifs[rbp]:
							out.write('{}\t{}\t{}\n'.format(motif, deg_motif, rbp))
						# print(motif, iter_items, deg_motif)
						# pass

	print("--- %s seconds ---" % (time.time() - start_time))
	print('DONE!', strftime('%a, %d %b %Y %I:%M:%S'))



if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Get sequences around SNPs')
	parser.add_argument('-i', metavar='indir', required=True, help='Input directory containing motifs')
	parser.add_argument('-r', metavar='rbns', required=True, help='RBNS 6mers')
	parser.add_argument('-o', metavar='outf', required=True, help='RBNS motifs found in our dataset')

	opts = parser.parse_args()
	print('Indir: %s' % opts.i)
	print('RBNS: %s' % opts.r)
	print('Outf: %s' % opts.o)

	main(sys.argv[1:])
