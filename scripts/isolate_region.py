import numpy as np
import argparse
import sys
try:
	from scripts.utils import import_aln
except:
	from utils import import_aln


def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--seqs', default='-')
	parser.add_argument('--ref', default=None)
	parser.add_argument('--region', nargs=2, type=int, default=None)
	args = parser.parse_args()
	#args.seqs = 'data/INSPIRE_071823_aln_gag_A1_bg_aln.fasta'
	#args.ref = "K03455.1"
	#args.region = [790,2292]
	if args.seqs == '-':
		seq_names, seqs = import_aln(sys.stdin)
	else:
		seq_names, seqs = import_aln(open(args.seqs, 'rt'))
	#seq_names, seqs = import_aln(args.seqs)
	seq_names_dict = {i:idx for idx, i in enumerate(seq_names)}
	if not args.ref:
		region = seqs[:,(args.region[0]-1):args.region[1]]
	if args.ref:
		ref_idx = [idx for idx, i in enumerate(seq_names) if args.ref in i]
		if len(ref_idx) > 1:
			raise Exception('reference found more than once in alignment')
		else:
			ref_idx = ref_idx[0]
		ref_cords = np.cumsum(seqs[ref_idx] != '-') -1
		region_start = \
			np.where(ref_cords == (args.region[0]))[0][0]
		region_end = \
			np.where(ref_cords == (args.region[1]))[0][0]
		region = seqs[:,(region_start-1):region_end]
	for idx, name in enumerate(seq_names):
		print('>'+name)
		print(''.join(region[idx]))



if __name__ == "__main__":
    run()

