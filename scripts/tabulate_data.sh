#!/usr/bin/env bash

#### -------------------------- ####
#### TABUALTE DESCRIPTIVE STATS ####
#### -------------------------- ####
# Demographics of study pop over time
# we pull ageyrs from internal file, not sharing publicly
Rscript scripts/tabulate_time_stratified_count.R \
	--dat data/rakai_drug_resistance_categorized_R15_R20_internal.tsv \
	--tableVars ageyrs \
	--out tmp1

Rscript scripts/tabulate_time_stratified_count.R \
	--dat data/rakai_drug_resistance_categorized_R15_R20.tsv \
	--tableVars age_cat comm_type sex \
	--out tmp2

cat tables/tmp1.tsv <(tail -n +3 tables/tmp2.tsv) > \
	tables/par_round_stratified_count.tsv

rm -rf tables/tmp*.tsv

# we pull ageyrs from internal file, not sharing publicly
Rscript scripts/tabulate_time_stratified_count.R \
	--dat data/rakai_drug_resistance_categorized_R15_R20_internal.tsv \
	--tableVars ageyrs \
	--out tmp1 \
	--filter "finalhiv == 'P'"

Rscript scripts/tabulate_time_stratified_count.R \
	--dat data/rakai_drug_resistance_categorized_R15_R20.tsv \
	--tableVars age_cat comm_type sex \
	--out tmp2 \
	--filter "finalhiv == 'P'"

cat tables/tmp1.tsv <(tail -n +3 tables/tmp2.tsv) > \
	tables/plhiv_round_stratified_count.tsv

rm -rf tables/tmp*.tsv

# Sequencing success among treatment-experienced viremic PLHIV 
# parsing of parenthesis here is weird
Rscript scripts/tabulate_multi_outcome.R \
	--dat data/rakai_drug_resistance_categorized_R15_R20.tsv \
	--filter 'finalhiv == "P" & viremic & round > 15 & !pre_treatment' \
	--tableVars round sequenced subtype_bestref_cat vl_copies_cat\
	--outcomes "!((is.na(sequence_id)))" \
	--out tmp1

Rscript scripts/tabulate_multi_outcome.R \
	--dat data/rakai_drug_resistance_categorized_R15_R20.tsv \
	--filter 'finalhiv == "P" & viremic & round > 15 & !pre_treatment & !is.na(sequence_id)' \
	--tableVars round sequenced subtype_bestref_cat  vl_copies_cat\
	--outcomes "!((is.na(insti)) | !(is.na(nnrti)) | !(is.na(nrti)) | !(is.na(pi))" "!(is.na(insti))" "!(is.na(nnrti))" "!(is.na(nrti))" "!(is.na(pi))" "!(is.na(insti)) & !(is.na(nnrti)) & !(is.na(nrti)) & !(is.na(pi)))"\
	--out tmp2

awk -F'\t' 'FNR==NR {vals[$1"_"$2] = $3"\t"$4; next}{print $1"\t"$2"\t"vals[$1"_"$2]"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' \
	tables/tmp1.tsv \
	tables/tmp2.tsv \
	> tables/sequenced_among_tx_viremic.tsv

rm -rf tables/tmp1.tsv
rm -rf tables/tmp2.tsv
# Sequencing success among pre-treatment viremic PLHIV 
# parsing of parenthesis here is weird
Rscript scripts/tabulate_multi_outcome.R \
	--dat data/rakai_drug_resistance_categorized_R15_R20.tsv \
	--filter 'finalhiv == "P" & viremic & pre_treatment' \
	--tableVars round sequenced subtype_bestref_cat vl_copies_cat\
	--outcomes "!((is.na(sequence_id)))" \
	--out tmp1

Rscript scripts/tabulate_multi_outcome.R \
	--dat data/rakai_drug_resistance_categorized_R15_R20.tsv \
	--filter 'finalhiv == "P" & viremic & pre_treatment & !is.na(sequence_id)' \
	--tableVars round sequenced subtype_bestref_cat  vl_copies_cat\
	--outcomes "!((is.na(insti)) | !(is.na(nnrti)) | !(is.na(nrti)) | !(is.na(pi))" "!(is.na(insti))" "!(is.na(nnrti))" "!(is.na(nrti))" "!(is.na(pi))" "!(is.na(insti)) & !(is.na(nnrti)) & !(is.na(nrti)) & !(is.na(pi)))"\
	--out tmp2

awk -F'\t' 'FNR==NR {vals[$1"_"$2] = $3"\t"$4; next}{print $1"\t"$2"\t"vals[$1"_"$2]"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' \
	tables/tmp1.tsv \
	tables/tmp2.tsv \
	> tables/sequenced_among_pt_viremic.tsv

rm -rf tables/tmp1.tsv
rm -rf tables/tmp2.tsv

#### ----------------------- ####
#### TABULATE MODEL OUTCOMES ####
#### ----------------------- ####
# Prevalence of HIV, tx HIV, and supp HIV by round
cat \
	<(head -n 1 models/hiv_among_par.tsv) \
	<(tail -n +2 models/hiv_among_par.tsv) \
	<(tail -n +2 models/tx_among_par.tsv) \
	<(tail -n +2 models/viremic_among_par.tsv) \
	> models/tmp_rr.tsv
cat \
	<(head -n 1 models/hiv_among_par_pred.tsv) \
	<(tail -n +2 models/hiv_among_par_pred.tsv) \
	<(tail -n +2 models/tx_among_par_pred.tsv) \
	<(tail -n +2 models/viremic_among_par_pred.tsv) \
	> models/tmp_prev_pred.tsv

Rscript scripts/tabulate_round_prev_rr.R \
	--prev models/tmp_prev_pred.tsv \
	--rr models/tmp_rr.tsv \
	--labelOrder hiv_among_par tx_among_par viremic_among_par \
	--nOut all \
	--out hiv_tx_viremic_among_par
rm -rf models/tmp_rr.tsv
rm -rf models/tmp_prev_pred.tsv

# Prevalence of suppression among PLHIV
cat \
	<(head -n 1 models/tx_among_plhiv.tsv) \
	<(tail -n +2 models/tx_among_plhiv.tsv) \
	<(tail -n +2 models/supp_among_plhiv.tsv) \
	> models/tmp_rr.tsv
cat \
	<(head -n 1 models/tx_among_plhiv_pred.tsv) \
	<(tail -n +2 models/tx_among_plhiv_pred.tsv) \
	<(tail -n +2 models/supp_among_plhiv_pred.tsv) \
	> models/tmp_prev_pred.tsv
Rscript scripts/tabulate_round_prev_rr.R \
	--prev models/tmp_prev_pred.tsv \
	--rr models/tmp_rr.tsv \
	--labelOrder tx_among_plhiv supp_among_plhiv \
	--nOut all \
	--out tx_supp_among_plhiv_round
rm -rf models/tmp_rr.tsv
rm -rf models/tmp_prev_pred.tsv

# Prevalence of suppression among tx
Rscript scripts/tabulate_round_prev_rr.R \
	--prev models/supp_among_tx_pred.tsv \
	--rr models/supp_among_tx.tsv \
	--labelOrder supp_among_tx \
	--nOut all \
	--out supp_among_tx_round


# Prevalence of resistance among treatment-experienced viremic PLHIV
# concatenate outputs
cat \
	models/any_among_viremic_tx_pred.tsv \
	<(tail -n +2 models/insti_among_viremic_tx_pred.tsv)\
	<(tail -n +2 models/nnrti_among_viremic_tx_pred.tsv)\
	<(tail -n +2 models/nrti_among_viremic_tx_pred.tsv)\
	<(tail -n +2 models/pi_among_viremic_tx_pred.tsv) \
	> models/tmp_pred.tsv

cat \
	models/any_among_viremic_tx.tsv \
	<(tail -n +2 models/insti_among_viremic_tx.tsv)\
	<(tail -n +2 models/nnrti_among_viremic_tx.tsv)\
	<(tail -n +2 models/nrti_among_viremic_tx.tsv)\
	<(tail -n +2 models/pi_among_viremic_tx.tsv) \
	> models/tmp_rr.tsv

Rscript scripts/tabulate_round_prev_rr.R \
	--prev models/tmp_pred.tsv \
	--rr models/tmp_rr.tsv \
	--labelOrder any_among_viremic_tx insti_among_viremic_tx nnrti_among_viremic_tx nrti_among_viremic_tx pi_among_viremic_tx\
	--out resistance_among_tx_round

# resistance risk factors among viremic tx
cat \
	<(awk -F'\t' '{OFS="\t"; $1=""; print "var"$0}' <(head -n 1 models/any_among_viremic_tx_sex_pred.tsv)) \
	<(tail -n +2 models/any_among_viremic_tx_sex_pred.tsv) \
	<(tail -n +2 models/any_among_viremic_tx_age_cat_pred.tsv) \
	<(tail -n +2 models/any_among_viremic_tx_comm_cat_pred.tsv) \
	<(tail -n +2 models/any_among_viremic_tx_subtype_bestref_cat_pred.tsv) \
	> models/tmp_pred.tsv

cat \
	models/any_among_viremic_tx_sex.tsv \
	<(tail -n +2 models/any_among_viremic_tx_age_cat.tsv) \
	<(tail -n +2 models/any_among_viremic_tx_comm_cat.tsv) \
	<(tail -n +2 models/any_among_viremic_tx_subtype_bestref_cat.tsv) \
	> models/tmp_rr.tsv

cat \
	models/any_among_viremic_tx_sex_round.tsv \
	<(tail -n +2 models/any_among_viremic_tx_age_cat_round.tsv) \
	<(tail -n +2 models/any_among_viremic_tx_comm_cat_round.tsv) \
	<(tail -n +2 models/any_among_viremic_tx_subtype_bestref_cat_round.tsv) \
	> models/tmp_rr_adj.tsv

Rscript scripts/tabulate_risk_factors.R \
	--out resistance_among_tx_covariates \
	--prev models/tmp_pred.tsv \
	--rr models/tmp_rr.tsv \
	--rrAdj models/tmp_rr_adj.tsv

# multi / single-class resistance
cat \
	models/multi-class_among_viremic_tx_pred.tsv \
	<(tail -n +2 models/single-class_among_viremic_tx_pred.tsv)\
	> models/tmp_pred.tsv

cat \
	models/multi-class_among_viremic_tx.tsv \
	<(tail -n +2 models/single-class_among_viremic_tx.tsv)\
	> models/tmp_rr.tsv

Rscript scripts/tabulate_round_prev_rr.R \
	--prev models/tmp_pred.tsv \
	--rr models/tmp_rr.tsv \
	--labelOrder multi-class_among_viremic_tx single-class_among_viremic_tx \
	--out multi_resistance_among_tx_round

# resistance mutations among viremic tx
cat \
	<(head -n 1 models/mutations/viremic-tx-inT97A_pred.tsv) \
	<(cat models/mutations/viremic-tx*_pred.tsv | grep -v "round\tfit") \
	> models/tmp_pred.tsv

cat \
	<(head -n 1 models/mutations/viremic-tx-inT97A.tsv) \
	<(cat models/mutations/viremic-tx*[^pred].tsv | grep -v "var\tRR") \
	> models/tmp_rr.tsv

Rscript scripts/tabulate_round_prev_rr.R \
	--prev models/tmp_pred.tsv \
	--rr models/tmp_rr.tsv \
	--nOut one \
	--labelOrder viremic-tx-inT97A viremic-tx-rtK103N viremic-tx-inS153Y viremic-tx-prL23I viremic-tx-rtM184V viremic-tx-rtM41L viremic-tx-rtP225H viremic-tx-rtT215Y viremic-tx-rtY181C viremic-tx-inE138K \
	--out mut_among_tx_round

# Prevalence of resistance among pre-treatment viremic PLHIV
# concatenate outputs
cat \
	models/insti_among_viremic_pt_pred.tsv \
	<(tail -n +2 models/nnrti_among_viremic_pt_pred.tsv)\
	<(tail -n +2 models/nrti_among_viremic_pt_pred.tsv)\
	<(tail -n +2 models/pi_among_viremic_pt_pred.tsv) \
	> models/tmp_pred.tsv

cat \
	models/insti_among_viremic_pt.tsv \
	<(tail -n +2 models/nnrti_among_viremic_pt.tsv)\
	<(tail -n +2 models/nrti_among_viremic_pt.tsv)\
	<(tail -n +2 models/pi_among_viremic_pt.tsv) \
	> models/tmp_rr.tsv

Rscript scripts/tabulate_round_prev_rr.R \
	--prev models/tmp_pred.tsv \
	--rr models/tmp_rr.tsv \
	--labelOrder insti_among_viremic_pt nnrti_among_viremic_pt nrti_among_viremic_pt pi_among_viremic_pt \
	--out resistance_among_pt_round

# resistance risk factors among viremic pt
cat \
	<(awk -F'\t' '{OFS="\t"; $1=""; print "var"$0}' <(head -n 1 models/nnrti_among_viremic_pt_sex_pred.tsv)) \
	<(tail -n +2 models/nnrti_among_viremic_pt_sex_pred.tsv) \
	<(tail -n +2 models/nnrti_among_viremic_pt_age_cat_pred.tsv) \
	<(tail -n +2 models/nnrti_among_viremic_pt_comm_cat_pred.tsv) \
	<(tail -n +2 models/nnrti_among_viremic_pt_subtype_bestref_cat_pred.tsv) \
	> models/tmp_pred.tsv

cat \
	models/nnrti_among_viremic_pt_sex.tsv \
	<(tail -n +2 models/nnrti_among_viremic_pt_age_cat.tsv) \
	<(tail -n +2 models/nnrti_among_viremic_pt_comm_cat.tsv) \
	<(tail -n +2 models/nnrti_among_viremic_pt_subtype_bestref_cat.tsv) \
	> models/tmp_rr.tsv

cat \
	models/nnrti_among_viremic_pt_sex_round.tsv \
	<(tail -n +2 models/nnrti_among_viremic_pt_age_cat_round.tsv) \
	<(tail -n +2 models/nnrti_among_viremic_pt_comm_cat_round.tsv) \
	<(tail -n +2 models/nnrti_among_viremic_pt_subtype_bestref_cat_round.tsv) \
	> models/tmp_rr_adj.tsv


Rscript scripts/tabulate_risk_factors.R \
	--out resistance_among_pt_covariates \
	--prev models/tmp_pred.tsv \
	--rr models/tmp_rr.tsv \
	--rrAdj models/tmp_rr_adj.tsv
	
# resistance mutations among viremic pt
cat \
	<(head -n 1 models/mutations/viremic-pt-rtK103N_pred.tsv) \
	<(cat models/mutations/viremic-pt*_pred.tsv | grep -v "round\tfit") \
	> models/tmp_pred.tsv

cat \
	<(head -n 1 models/mutations/viremic-pt-rtK103N.tsv) \
	<(cat models/mutations/viremic-pt*[^pred].tsv | grep -v "var\tRR") \
	> models/tmp_rr.tsv

Rscript scripts/tabulate_round_prev_rr.R \
	--prev models/tmp_pred.tsv \
	--rr models/tmp_rr.tsv \
	--nOut one \
	--labelOrder viremic-pt-rtK103N viremic-pt-inT97A viremic-pt-inS153Y viremic-pt-rtE138A viremic-pt-prL23I viremic-pt-prL33F viremic-pt-prG73V viremic-pt-rtF77L viremic-pt-rtG190A viremic-pt-rtM41L \
	--out mut_among_pt_round

# finally, insti resistance
Rscript scripts/tabulate_insti_resistance.R