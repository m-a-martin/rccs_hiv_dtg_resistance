suppressMessages(library(tidyverse))
suppressMessages(library(geepack))
suppressMessages(library(emmeans))
suppressMessages(library(patchwork))
source('scripts/utils.R')
 

# todo move all of this data processing to data processing script
d = read_tsv('data/rakai_drug_resistance_categorized_R15_R20.tsv', show_col_types=FALSE)  %>% 
	mutate(
		numeric_round = as.numeric(round),
		round = as.character(round),
		round_year = str_split(round_mid_date, "-", simplify=TRUE)[,1])

hiv_dr_mut = read_tsv('data/rakai_drug_resistance_mut_R15_R20.tsv',
	show_col_types=FALSE) %>%
	filter(iteration == '5pct_10' & mut != "WT") %>%
	mutate(round = as.character(round))

# those participants with intermediate/high INSTI resistance
dtg = d %>% filter(DTG == '2' | DTG == '3' | DTG == '4') %>%
	select(study_id, sex, age_cat, round, round_year, Mut, pre_treatment, first_arv, first_arv_date, first_pos_date, last_neg_date) %>%
	arrange(round, -pre_treatment)


# add first arv dates
#dtg = dtg %>%
#	left_join(
#		read_csv('data/R20_15_consolidated.csv', show_col_types=FALSE) %>% 
#			select(study_id, first_arv, first_arv_date) %>% unique(),
#		by=c('study_id', 'first_arv'))

# format muts
muts = sapply(str_split(dtg$Mut, '\\|'), function(x){str_split(x, ";", simplify=TRUE)[,1]})
in_muts = sapply(muts, function(x){paste(x[grepl('in', x)], collapse=', ')})
other_muts = sapply(muts, function(x){paste(x[!grepl('in', x)], collapse=', ')})
dtg$insti_mutations = in_muts
dtg$other_mutations = other_muts

# add treatment status
dtg = dtg %>% mutate(
	treatment_status = if_else(pre_treatment == TRUE, 'pre-treatment', 'treatment-experienced'),
	first_arv = if_else(first_arv == Inf, NA, first_arv))

# converts dates to month years
convert_date_fmt = function(x){
	format(
		as.Date(
			if_else(is.na(x), NA, paste0(x, '-01')),
			format='%Y-%m-%d'),
		'%B %Y')
}


dtg = dtg %>% mutate(
	first_arv_date = convert_date_fmt(first_arv_date),
	last_neg_date = convert_date_fmt(last_neg_date),
	first_pos_date = convert_date_fmt(first_pos_date)) %>% print(width=Inf)

# choose columns
# sex	age	round	INSTI mutations	other DR mutations	pre_treatment	last HIV- round first HIV+ round	
write_tsv(
	dtg %>% select(sex, age_cat, round_year, insti_mutations, other_mutations, treatment_status, 
		first_arv_date, last_neg_date, first_pos_date),
	'tables/insti_mutations.tsv',
	na="")

write_tsv(
	dtg %>% select(study_id, sex, age_cat, round_year, insti_mutations, other_mutations, treatment_status, 
		first_arv_date, last_neg_date, first_pos_date),
	'tables/insti_mutations_internal.tsv',
	na="")
