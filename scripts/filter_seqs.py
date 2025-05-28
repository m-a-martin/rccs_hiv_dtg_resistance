import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import sys
try:
	from scripts.utils import import_aln, plot_style
except:
	from utils import import_aln, plot_style



def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--seqs', default=None)
	parser.add_argument('--ref', default=None)
	parser.add_argument('--minACTG', default=1, type=float)
	args = parser.parse_args()
	#args.seqs = 'data/bg_A1_pol_aln.fasta'
	#args.minACTG = 0.95
	#args.ref = "K03455.1"

	if args.seqs == '-':
		seq_names, seqs = import_aln(sys.stdin)
	else:
		seq_names, seqs = import_aln(open(args.seqs, 'rt'))

	ref_idx = [idx for idx, i in enumerate(seq_names) if args.ref in i]
	if len(ref_idx) > 1:
		raise Exception('reference found more than once in alignment')
	else:
		ref_idx = ref_idx[0]
	
	p_valid = np.isin(seqs[:,seqs[ref_idx] != '-'], ['A', 'a', 'C', 'c', 'T', 't', 'G', 'g']).sum(axis=1)/(seqs[ref_idx] != '-').sum()
	for n,s,p in zip(seq_names, seqs, p_valid):
		if p > args.minACTG:
			print('>'+n)
			print(''.join(s))


if __name__ == "__main__":
    run()










