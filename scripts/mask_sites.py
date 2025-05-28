import numpy as np
import pandas as pd
from collections import Counter
import argparse
try:
	from scripts.utils import import_aln
except:
	from utils import import_aln


from collections import defaultdict
codon_dict = defaultdict(lambda: np.nan)
codon_dict.update(
	{'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 
	'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 
	'TAT': 'Y', 'TAC': 'Y', 'TAA': 'X', 'TAG': 'X', 
	'TGT': 'C', 'TGC': 'C', 'TGA': 'X', 'TGG': 'W', 
	'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 
	'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 
	'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 
	'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 
	'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 
	'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 
	'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 
	'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 
	'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 
	'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 
	'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 
	'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'})


def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--seqs', default=None)
	parser.add_argument('--maskSites', default=None)
	parser.add_argument('--maskClassCol', default=0, type=int)
	parser.add_argument('--maskSiteCol', default=1, type=int)
	parser.add_argument('--genePos', default=None)
	parser.add_argument('--genePosClassCol', default=0, type=int)
	parser.add_argument('--genePosStartCol', default=2, type=int)
	parser.add_argument('--ref', type=str)
	args = parser.parse_args()
	#args.seqs = 'data/bg_A1_pol_realn.fasta'
	#args.genePos = 'data/class_gene_position.csv'
	#args.maskSites = 'data/dr_mutations.csv'
	#args.ref = "K03455.1"
	seq_names, seqs = import_aln(open(args.seqs, 'r'))
	ref_idx = [idx for idx, i in enumerate(seq_names) if args.ref in i]
	if len(ref_idx) > 1:
		raise Exception('reference found more than once in alignment')
	else:
		ref_idx = ref_idx[0]
	gene_pos = pd.read_csv(args.genePos)
	mask_sites = pd.read_csv(args.maskSites, sep=',', low_memory=False)
	mask_sites = mask_sites.merge(gene_pos, on='class').\
		assign(nuc_pos = lambda k: k.site*3 + (k.start-1))
	dg_seqs = seqs[:,seqs[ref_idx] != '-']
	keep = np.array([True]*dg_seqs.shape[1])
	for idx, row in mask_sites.iterrows():
		# first confirm that wt AA is the most common
		# get codon
		codon_cords = np.arange(row.nuc_pos-3, row.nuc_pos)
		codons = Counter([''.join(i) for i in dg_seqs[:,codon_cords]])
		inferred_wt_aa = codon_dict[list(codons.keys())[0].upper()]
		if inferred_wt_aa != row.wt:
			raise Exception('wild-type amino acids do not match, are you sure your numbering is correct?')
		keep[codon_cords] = False

	masked_dg_seqs = dg_seqs[:,keep]
	for name, seq in zip(seq_names, masked_dg_seqs):
		print('>' + name)
		print(''.join(seq))

	

if __name__ == "__main__":
    run()

