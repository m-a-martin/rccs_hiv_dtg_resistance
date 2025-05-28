import argparse
import numpy as np
import pandas as pd
from itertools import product
try:
	from scripts.utils import import_aln, map_arr
except:
	from utils import import_aln, map_arr


def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--seqs', default=None)
	parser.add_argument('--dists', default=None)
	parser.add_argument('--dats', default=None, nargs='+')
	parser.add_argument('--muts', default=None, nargs='+')
	args = parser.parse_args()
	#args.seqs = 'data/2024-10-02_pangea2_aln_pol_mask.fasta'
	#args.dists = 'data/2024-10-02_pangea2_aln_pol_mask_dist.csv'
	#args.dats = ['data/rakai_drug_resistance_categorized_R15_R20.tsv', 'data/other_rakai_drug_resistance_categorized.tsv']
	#args.muts = ['data/rakai_drug_resistance_mut_R15_R20.tsv', 'data/other_rakai_drug_resistance_mut.tsv']
	dat = pd.concat([
		pd.read_csv(i, sep='\t')[['study_id', 'int_date', 'subtype_bestref']] for i in args.dats]).\
		assign(label = lambda k: k.study_id.astype(str) + '_' + \
			k.int_date.str.split('T', expand=True)[0])

	seq_names, seqs = import_aln(open(args.seqs, 'r'))
	seq_names_dict = {i:idx for idx, i in enumerate(seq_names)}
	dists = pd.read_csv(args.dists).\
		query('seq1 != seq2').\
		query("seq1 != 'K03455.1' & seq2 != 'K03455.1'")
	# set diag to Inf
	muts = pd.concat([
		pd.read_csv(i, sep='\t') for i in args.muts])
	# get inS153Y seqs
	inS153Y_seqs = muts.query('mut == "inS153Y"').\
		assign(label = lambda k: k['study_id'].astype(str) + '_' + \
			k['int_date'].astype(str).str.split('T', expand=True)[0]).label.drop_duplicates().values

	# filter dists for only intra-subtype
	# and sampled post-R20
	st_dists = dists.\
		merge(
			dat[['label', 'subtype_bestref', 'int_date']].\
			rename(columns={
				'label': 'seq2',
				'subtype_bestref': 'subtype_bestref2',
				'int_date': 'int_date2'}),
			how='left').\
		merge(
			dat[['label', 'subtype_bestref', 'int_date']].\
			rename(columns={
				'label': 'seq1',
				'subtype_bestref': 'subtype_bestref1',
				'int_date': 'int_date1'}),
			how='left').\
		query('subtype_bestref1 == subtype_bestref2')
	##& int_date1 > "2021-02-01" & int_date2 > "2021-02-01"
	# had this code to filter by subtype
	# but turns out there just aren't that many r20 sequences
	# so including all! 
	use_seqs = st_dists.loc[np.isin(st_dists.seq1, inS153Y_seqs),:].\
		sort_values(by='d').\
		groupby('seq1').\
		head(10).\
		query('subtype_bestref1 == "A1" | subtype_bestref1 == "D"').\
		drop(['subtype_bestref2', 'd', 'int_date1', 'int_date2'], axis=1).\
		melt(id_vars='subtype_bestref1').\
		drop(['variable'], axis=1).\
		drop_duplicates().\
		assign(
			seq = lambda k: [''.join(i) for i in seqs[k.value.map(seq_names_dict)]],
			value = lambda k: k.value + '_' + np.isin(k.value, inS153Y_seqs).astype(str))

	# splt by subtype
	for g, g_dat in use_seqs.groupby('subtype_bestref1'):
		with open('.'.join(args.seqs.split('.')[:-1]) + '_' + g + '_inS153Y.fasta', 'w') as fp:
			fp.write('>' + seq_names[0].split(' ')[0] + '\n')
			fp.write(''.join(seqs[0,:]) + '\n')
			for idx, i in g_dat.iterrows():
				fp.write('>' + i.value + '\n')
				fp.write(i.seq + '\n')



if __name__ == "__main__":
    run()

