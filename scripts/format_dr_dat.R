suppressMessages(require(tidyverse))
suppressMessages(require(readxl))
suppressMessages(require(haven))
suppressMessages(source('scripts/utils.R'))

get_resistance = function(.data, cols, label){
	.data = .data %>% mutate(!!label := case_when(
		# if any of the columns are a 3 and 4 and 
		# all of the columns are in S, 1, 2, 3, or 4 (i.e. not X) then call it intermediate/high
		if_any(all_of(cols), function(x) x %in% c("3", "4")) & 
			if_all(all_of(cols), function(x) x %in% c("P", "S", "0", "1", "2", "3", "4")) ~ 'intermediate/high',
		# if all of the columns are an 'S', '1', or '2'
		# then call it susceptible
		if_all(all_of(cols), function(x) x %in% c("P", "S", "0", "1", "2")) ~ 'susceptible'))
		# otherwise will return NA
	return(.data)
}

rccs_metadata = read_csv('data/R20_15_consolidated.csv', show_col_types=FALSE)

# create dummy study id variable
# need to get all study IDs that are used across RCCS rounds 15-20 and non-RCCS sequenced samples
study_id_map = bind_rows(
		rccs_metadata %>%
			select(study_id) %>%
			unique(),
		read_csv('data/2024-10-02_pangea2_mike_sequence_id_metadata.csv', show_col_types=FALSE) %>%
			# remove "RK" so matches raka metadata
			mutate(
				study_id = str_split(pt_id, "-", simplify=TRUE)[,2]) %>%
			select(study_id) %>%
			unique()) %>%
	unique() %>%
	mutate(
		study_id_internal = study_id,
		study_id = sprintf(paste0('%0',  nchar(as.character(n())), 'd'), seq(1,n())))

write_tsv(study_id_map, 'data/study_id_map.tsv')

rccs_metadata = rccs_metadata %>% 
	rename('study_id_internal'='study_id') %>%
	left_join(study_id_map, 'study_id_internal')

# get survey dates through R19
survey_dates = rccs_metadata %>% group_by(round) %>% 
	summarise(
		round_min_date = quantile(int_date, c(0.01), type=1, na.rm=TRUE),
		round_mid_date = quantile(int_date, c(0.5), type=1, na.rm=TRUE),
		round_max_date = quantile(int_date, c(0.99), type=1, na.rm=TRUE)
		)

# get sequence data
format_seq_dat = read_csv('data/2024-10-02_pangea2_mike_sequence_id_metadata.csv', show_col_types=FALSE) %>%
	# remove "RK" so matches raka metadata
	mutate(
		study_id = str_split(pt_id, "-", simplify=TRUE)[,2],
		sequenced = if_else(substring(sequence_id, 1, 2) == 'PG', 'PANGEA1', 'PANGEA2')) %>%
	# rename so matches rakai data 
	rename(
		sequence_id_internal=sequence_id,
		int_date=visit_dt,
		study_id_internal=study_id) %>%
	# keeping length so we can filter on sequencing quality
	# for samples that were resequenced
	select(sequence_id_internal, study_id_internal, sex, int_date, length_strict, readnum_hiv, subtype_bestref, shiver_consensus, sequenced, vl_dt, vl_result) %>%
	unique() %>%
	left_join(study_id_map, by='study_id_internal') %>%
	rowwise() %>%
	mutate(sequence_id = sub(study_id_internal, study_id, sequence_id_internal))

write_tsv(format_seq_dat %>% relocate(study_id) %>% relocate(sequence_id), 
	'data/2024-10-02_pangea2_mike_sequence_id_metadata_internal.tsv')

id_map = format_seq_dat %>% select(-shiver_consensus)

hiv_dr = 
	bind_rows(
		read_csv("data/HIVdbOUT_DRMpredictions_Summary_2024-06-13_drmSEQ_Rakai.csv", show_col_types=FALSE),
		read_csv("data/HIVdbOUT_DRMpredictions_Summary_2024-09-09_drmSEQ_Rakai_failed_samples.csv", show_col_types=FALSE)) %>%
	mutate(
		# cut-off used to call resistance
		iteration=str_replace(iteration, "remap_", "")) %>%
	rename(sequence_id_internal = sequence_id)

names(hiv_dr) = gsub('/r', '', names(hiv_dr))

hiv_dr$n_X = rowSums(hiv_dr[,c(
	insti_cols,
	nnrti_cols,
	nrti_cols, 
	pi_cols)] == 'X')

hiv_dr = hiv_dr %>% 
	filter(iteration == "5pct_10") %>%	
	left_join(id_map, by=c('sequence_id_internal'), relationship='many-to-many') %>%
	group_by(study_id_internal, int_date) %>%
	filter(readnum_hiv == max(readnum_hiv)) %>%
	filter(n_X == min(n_X)) %>%
	slice(1) 

hiv_dr = hiv_dr %>% 
			mutate(
				dr_dat = TRUE,
				valid_dr_dat = !if_all(all_of(all_cols), function(x) x %in% c("X")))

# categorize
# todo make this faster
# keep DTG column
hiv_dr_cat = hiv_dr %>%
	get_resistance(c('DTG'), 'dtg') %>% 
	get_resistance(insti_cols, "insti") %>%
	get_resistance(nnrti_cols, "nnrti") %>%
	get_resistance(nrti_cols, "nrti") %>%
	get_resistance(pi_cols, "pi") %>%
	select(-all_of(c(insti_cols[!(insti_cols == 'DTG')], nnrti_cols, nrti_cols, pi_cols))) %>%
	select(-c(Haplo, n_X))

# merge with rccs
rccs_hiv_dr_cat = hiv_dr_cat %>%
	select(-sex) %>%
	right_join(rccs_metadata, by=c('study_id', 'study_id_internal', 'int_date')) %>%
	mutate(
		dr_dat = case_when(
			finalhiv == 'P' & !is.na(dr_dat) ~ dr_dat,
			finalhiv == 'P' & is.na(dr_dat) ~ FALSE,
			TRUE ~ dr_dat),
		valid_dr_dat = case_when(
			finalhiv == 'P' & !is.na(valid_dr_dat) ~ valid_dr_dat,
			finalhiv == 'P' & is.na(valid_dr_dat) ~ FALSE,
			TRUE ~ valid_dr_dat))

# impute R15 pre-treatment viral loads
# among pre-treatment PLHIV only
# include those who are self-reported at pre-treatment but have suppressed VL in R15
set.seed(1)
prob_viremic = rccs_hiv_dr_cat %>% 
	filter(
		round == 15 & 
		finalhiv == "P" & 
		first_arv > 15 &
		first_supr >= 15 &
		!is.na(numeric_copies) & 
		!is.na(valid_dr_dat)) %>%
	group_by(round, valid_dr_dat) %>%
	summarise(p_viremic = sum(numeric_copies > 1000)/n(), .groups="drop") %>%
	ungroup() %>%
	mutate(pre_treatment = TRUE)

rccs_hiv_dr_cat = rccs_hiv_dr_cat %>% 
	left_join(
		prob_viremic,
		by=c('round', 'pre_treatment', 'valid_dr_dat')) %>%
	mutate(
		viremic_raw = numeric_copies > 1000,
		viremic = case_when(
			numeric_copies > 1000 & finalhiv == 'P' ~ TRUE,
			numeric_copies <= 1000 & finalhiv == 'P' ~ FALSE,
			(round == 15) & pre_treatment == TRUE & is.na(numeric_copies) & finalhiv == 'P' ~  
				as.logical(rbinom(n(), 1, replace_na(p_viremic, 0)))))

# categorize subtype
rccs_hiv_dr_cat = rccs_hiv_dr_cat %>%
	ungroup() %>% 
	mutate(subtype_bestref_cat = case_when(
		is.na(subtype_bestref) ~ NA,
		subtype_bestref == 'A1' ~ subtype_bestref,
		subtype_bestref == 'D' ~ subtype_bestref,
		subtype_bestref == 'C' ~ subtype_bestref,
		TRUE ~ 'other'))

write_tsv(rccs_hiv_dr_cat %>% select(-vl_dt, -vl_result, -study_id_internal, -sequence_id_internal) %>%
		relocate(study_id) %>%
		relocate(sequence_id), 
	'data/rakai_drug_resistance_categorized_R15_R20.tsv')

write_tsv(rccs_hiv_dr_cat %>% select(-vl_dt, -vl_result), 
	'data/rakai_drug_resistance_categorized_R15_R20_internal.tsv')


# get non cohort categorized data
# filter by minimum cohort date
# drop those with viral loads from earlier rounds
rakai_hiv_dr_cat = hiv_dr_cat %>%
	ungroup() %>%
	left_join(rccs_metadata %>% select(study_id_internal, int_date) %>% mutate(cohort = TRUE), 
		by=c('study_id_internal', 'int_date')) %>%
	filter(is.na(cohort)) %>%
	select(-cohort) %>%
	filter(int_date >= (survey_dates %>% filter(round == 15))$round_min_date & 
		int_date == vl_dt) %>%
	rename(numeric_copies = vl_result) %>%
	mutate(viremic = numeric_copies > 1000) %>%
	left_join(study_id_map,
		by=c('study_id', 'study_id_internal'))

write_tsv(rakai_hiv_dr_cat %>% 
		select(-study_id_internal, -sequence_id_internal) %>%
		relocate(study_id) %>%
		relocate(sequence_id),
		'data/other_rakai_drug_resistance_categorized.tsv')

write_tsv(rakai_hiv_dr_cat,
	'data/other_rakai_drug_resistance_categorized_internal.tsv')

# now format the haplotype column
# add number of mutations
hiv_dr = hiv_dr %>% mutate(n_muts = replace_na(str_count(Mut, "\\|")+1, 0))

extract = function(x, idx){
	if (all(is.na(x))){
		return(NA)
	}else{
		return(x[idx])
	}
}


extract_col = function(x, idx){
	if (all(is.na(x))){
		return(NA)
	}else{
		return(x[,idx])
	}
}


hiv_muts = bind_cols(
	tibble(
		mut = 
			replace_na(
				unlist(lapply(
					str_split(hiv_dr$Mut, "\\|"),
					function(x) str_split(x, ";", simplify=TRUE)[,1])),
				'WT'),
		freq = 
			replace_na(
				unlist(
					lapply(str_split(hiv_dr$Mut, "\\|"), function(x) extract_col(str_split(x, ";", simplify=TRUE), 3))),
				""),
		reads = 
			replace_na(
				unlist(
					lapply(str_split(hiv_dr$Mut, "\\|"), function(x) extract_col(str_split(x, ";", simplify=TRUE), 2))),
				"")),	
	(hiv_dr %>% 
				select(study_id_internal, study_id, sequence_id_internal, sequence_id, int_date, iteration))[
			rep(seq(1,nrow(hiv_dr)), if_else(hiv_dr$n_muts==0, 1, hiv_dr$n_muts)),])


write_tsv(
	hiv_muts %>%
		inner_join(rccs_metadata %>% select(study_id_internal, int_date, round),
			by=c('study_id_internal', 'int_date')) %>%
		select(-study_id_internal, -sequence_id_internal), 
	'data/rakai_drug_resistance_mut_R15_R20.tsv')

write_tsv(
	hiv_muts %>%
		inner_join(rccs_metadata %>% select(study_id_internal, int_date, round),
			by=c('study_id_internal', 'int_date')), 
	'data/rakai_drug_resistance_mut_R15_R20_internal.tsv')


write_tsv(
	hiv_muts %>%
		#left_join(study_id_map,
		#		by=c('study_id_internal')) %>%
			left_join(rccs_metadata %>% select(study_id_internal, study_id, int_date, round),
				by=c('study_id_internal', 'study_id', 'int_date')) %>%
			filter(is.na(round)) %>%
			select(-round, -study_id_internal, -sequence_id_internal), 
	'data/other_rakai_drug_resistance_mut.tsv')


write_tsv(
	hiv_muts %>%
		#left_join(study_id_map,
		#		by=c('study_id_internal')) %>%
			left_join(rccs_metadata %>% select(study_id_internal, study_id, int_date, round),
				by=c('study_id_internal', 'study_id', 'int_date')) %>%
			filter(is.na(round)) %>%
			select(-round), 
	'data/other_rakai_drug_resistance_mut_internal.tsv')


#rccs_metadata = read_csv('data/R20_15_consolidated.csv', show_col_types=FALSE)
#old_rccs_metadata = read_format_rccs_dat('data/RCCSdata_R001_R019_VOIs.dta') 
#all_rccs_metadata = bind_rows(
#	rccs_metadata %>% select(study_id, int_date, round) %>%
#		mutate(int_date = as.Date(int_date)),
#	old_rccs_metadata %>% filter(round < 15) %>% select(study_id, int_date, round))

#write_tsv(  
#	hiv_muts %>%
#		left_join(all_rccs_metadata %>% select(study_id, int_date, round), by=c('study_id', 'int_date')), 
#	'data/rakai_drug_resistance_mut_all.tsv')

