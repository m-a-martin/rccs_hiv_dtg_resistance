#!/usr/bin/env bash

# filter metadata file for those just in categorized resistance output
awk -F'\t' 'FNR==NR{a[$1] = "true";next}{if ($1 in a) print $0}'\
	<(cat \
		data/rakai_drug_resistance_categorized_R15_R20.tsv \
		data/other_rakai_drug_resistance_categorized.tsv) \
	data/2024-10-02_pangea2_mike_sequence_id_metadata_internal.tsv \
	> data/2024-10-02_pangea2_mike_sequence_id_metadata_res.tsv

# convert metadata file to fasta file
awk -F'\t' 'FNR==NR{for(i=1; i<=NF; i++) {cols[$i] = i}; next}{print ">"$cols["study_id"]"_"$cols["int_date"]"\n"$cols["shiver_consensus"]}' \
	<(head -n 1 data/2024-10-02_pangea2_mike_sequence_id_metadata_res.tsv) \
	<(tail -n +2 data/2024-10-02_pangea2_mike_sequence_id_metadata_res.tsv) \
	> data/2024-10-02_pangea2.fasta

# align seqs
mafft --thread 4 --auto \
	<(cat \
			data/K03455.1.fasta \
			data/2024-10-02_pangea2.fasta) > data/2024-10-02_pangea2_aln.fasta

# isolate polymerase and filter for completeness
python3 scripts/isolate_region.py \
	--seqs data/2024-10-02_pangea2_aln.fasta \
	--ref "K03455.1" \
	--region 2085 5096 | \
python3 scripts/filter_seqs.py \
	--seqs - \
	--ref "K03455.1" \
	--minACTG 0.95 \
> data/2024-10-02_pangea2_aln_pol.fasta

# mask resistance sites
python3 scripts/mask_sites.py \
	--seqs data/2024-10-02_pangea2_aln_pol.fasta \
	--maskSites data/dr_mutations.csv \
	--genePos data/class_gene_position.csv \
	--ref "K03455.1" \
	> data/2024-10-02_pangea2_aln_pol_mask.fasta

# calculate pairwise distnaces
python3 scripts/calc_distance.py \
	--seqs data/2024-10-02_pangea2_aln_pol_mask.fasta

# get tree seqs
cat \
	'data/rakai_drug_resistance_categorized_R15_R20.tsv' \
	<(tail -n +2 'data/other_rakai_drug_resistance_categorized.tsv') \
	> data/tmp.tsv

cat \
	'data/rakai_drug_resistance_mut_R15_R20.tsv' \
	<(tail -n +2 'data/other_rakai_drug_resistance_mut.tsv') \
	> data/tmp2.tsv

python3 scripts/get_tree_seqs.py \
	--seqs data/2024-10-02_pangea2_aln_pol_mask.fasta \
	--dists data/2024-10-02_pangea2_aln_pol_mask_dist.csv \
	--dats data/rakai_drug_resistance_categorized_R15_R20.tsv data/other_rakai_drug_resistance_categorized.tsv \
	--muts data/rakai_drug_resistance_mut_R15_R20.tsv data/other_rakai_drug_resistance_mut.tsv

iqtree2 -s data/2024-10-02_pangea2_aln_pol_mask_A1_inS153Y.fasta \
	-o K03455.1
iqtree2 -s data/2024-10-02_pangea2_aln_pol_mask_D_inS153Y.fasta \
	-o K03455.1
