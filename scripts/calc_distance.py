import argparse
import numpy as np
from itertools import product
try:
	from scripts.utils import import_aln, map_arr
except:
	from utils import import_aln, map_arr


def k2p_d(s1,s2, gt1, gt2, n1, n2):
	# how many sites are genotyped in both
	n = gt1 & gt2
	n_sum = n.sum()
	# which sites have differences
	d = np.where((s1 != s2) & n)
	# for each difference
	# is this difference a transition (FALSE) 
	# or transversion (TRUE)
	d_type = n1[d[0]] == n2[d[0]]
	p = d_type.sum()/n_sum
	q = (~d_type).sum()/n_sum
	return(-(1/2)*np.log((1-2*p-q)*np.sqrt(1-2*q)))


def k2p_arr_d(seq_arr):
	# https://doi.org/10.1007/BF01731581
	# d = -(1/2)ln((1-2P-Q)sqrt(1-2Q))
	# p = transition frequency
	# q = transversion frequency
	# dictionary that maps nucleotides to purine (0) or pyrimidine (1) state
	from collections import defaultdict
	n_type_dict = defaultdict(lambda: np.nan)
	n_type_dict.update({
		'a': 0,
		'c': 1,
		'g': 0,
		't': 1})
	# clean input array to just unambiguously genotyped nucleotides 
	# this collapses gaps ("-") to n, which are distinct states, but 
	# both are exlcuded for distance calculation, so fine here
	valid = np.array(['a', 'c', 't', 'g'])
	gt = np.isin(seq_arr, ['a', 'c', 't', 'g'])
	seq_arr_clean = np.where(~gt, 'n', seq_arr)
	# for each item in seq_arr is it purine or pyridime 
	n_type = map_arr(seq_arr_clean, n_type_dict)
	# empty distance arr
	d_arr = np.zeros((seq_arr.shape[0], seq_arr.shape[0]))
	# iterate over sequences in the input array and compare them to every other sequence
	p = np.zeros(seq_arr.shape[0])
	q = np.zeros(seq_arr.shape[0])
	for i in np.arange(seq_arr.shape[0]):
		d_arr[i,i] = 0
		for j in np.arange(i+1,seq_arr.shape[0]):
			i_j_d = k2p_d(seq_arr_clean[i],
				seq_arr_clean[j],
				gt[i],
				gt[j],
				n_type[i],
				n_type[j])
			d_arr[i,j] = i_j_d
			d_arr[j,i] = i_j_d
	return(d_arr)



def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--seqs', default=None)
	args = parser.parse_args()
	#args.seqs = 'data/2024-10-02_pangea2_aln_pol_mask.fasta'
	out_name = '.'.join(args.seqs.split('.')[:-1]) + '_dist.csv'
	seq_names, seqs = import_aln(open(args.seqs, 'r'))

	#seqs = seqs[:100,:]	
	#seq_names = seq_names[:10]
	d_arr = k2p_arr_d(seqs)


	with open(out_name, 'w') as fp:
		fp.write('seq1,seq2,d\n')
		for n, d in zip(product(seq_names, repeat=2), d_arr.flatten()):
			_ = fp.write(f'{n[0].split(" ")[0]},{n[1].split(" ")[0]},{d}\n')



if __name__ == "__main__":
    run()


