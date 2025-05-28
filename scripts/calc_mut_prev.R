suppressMessages(library(tidyverse))
suppressMessages(library(geepack))
suppressMessages(library(emmeans))
suppressMessages(library(patchwork))
suppressMessages(source('scripts/utils.R'))
 

# todo move all of this data processing to data processing script
d = read_tsv('data/rakai_drug_resistance_categorized_R15_R20.tsv', show_col_types=FALSE)  %>%
	group_by(round) %>%
	mutate(
		median_int_date = median(int_date, na.rm=TRUE)) %>%
	ungroup() %>%
	mutate(
		numeric_round = as.numeric(round),
		round = as.character(round),
		pre_treatment = round < first_arv,
		numeric_copies = case_when(
			is.na(copies) ~ NA,
			copies == 'INV.IC' ~ NA,
			copies == 'BD' ~ -1,
			copies == 'ND' ~ -1,
			copies == 'QNS' ~ -1,
			grepl('<', copies) ~ -1,
			grepl('Not Dete', copies) ~ -1,
			TRUE ~ as.numeric(copies)))


hiv_dr_mut = read_tsv('data/rakai_drug_resistance_mut_R15_R20.tsv',
	show_col_types=FALSE) %>%
	filter(iteration == '5pct_10' & mut != "WT") %>%
	mutate(round = as.character(round))


filters = c(
	"(round > 15) & (finalhiv == 'P') & (numeric_copies > 1000) & (pre_treatment == FALSE)",
	"(finalhiv == 'P') & (numeric_copies > 1000) & (pre_treatment == TRUE)")

labels = c(
	'viremic-tx',
	'viremic-pt')

min_outcomes = Inf
for (i in 1:length(filters)){
	filter = filters[i]
	i_d = d %>% filter(eval(parse(text=filter)))
	muts = i_d %>% select(study_id, round, median_int_date) %>%
		inner_join(hiv_dr_mut, by=c('study_id', 'round'))
	get_muts = muts %>% filter(round == max(as.numeric(round))) %>%
		group_by(mut) %>%
		summarise(n=n()) %>%
		arrange(-n) %>%
		slice(1:10)
	for (idx in 1:nrow(get_muts)){
		i_mut = get_muts$mut[idx]
		print(i_mut)
		out_file = paste(c('models/mutations/', labels[i], '-', i_mut), collapse='')
		m_d = i_d %>%
			left_join(
				hiv_dr_mut %>% filter(mut ==i_mut) %>%
					select(mut, study_id, round) %>%
					mutate(y = TRUE),
					by=c('study_id', 'round')) %>%
			mutate(y = case_when(
				!is.na(y) & (y == TRUE) & (!is.na(insti) & !is.na(nnrti) & !is.na(nrti) & !is.na(pi)) ~ TRUE,
				!is.na(y) & (y == TRUE) & (is.na(insti) | is.na(nnrti) | is.na(nrti) | is.na(pi)) ~ NA,
				is.na(y) & (!is.na(insti) & !is.na(nnrti) & !is.na(nrti) & !is.na(pi)) ~ FALSE,
				is.na(y) & (is.na(insti) | is.na(nnrti) | is.na(nrti) | is.na(pi)) ~ NA))
		# fit model
		m_output = run_process_model(
			m_d, 
			as.formula(if_else(length(unique(m_d$round)) > 1, 'y ~ round', 'y ~ 1')), 
			out_file,
			corstr=if_else(sum(m_d$y, na.rm=TRUE) > min_outcomes, 'fit', 'independence'), 
			id_col='study_id',
			order_col='round')
		# output model summary
		m_summary = as_tibble(m_output$o, rownames='var') %>% 
				mutate(
					label = paste(labels[i], i_mut, sep='-'),
					corr = m_output$min_qic_corr)
		# add median int dates if round in pred
		if ('round' %in% colnames(m_summary)){
			m_summary = m_summary %>% left_join(m_d %>% select(round, median_int_date) %>% unique(), by='round')
		}
		write_tsv(m_summary, 
			paste(out_file, '.tsv', sep=''))
		# predict values in each round
		pred = as_tibble(emmeans(
					m_output$m,  
					as.formula(if_else(length(unique(m_d$round)) > 1, ' ~ round', ' ~ 1')),
					type="response")) %>%
					rename(c(
						'fit' = 'rate',
						'se.fit' = 'SE',
						'lwr'= 'lower.CL', 
						'upr'='upper.CL'))
		if (length(unique(m_d$round)) > 1){
			pred = pred %>% left_join(m_d %>% group_by(round, median_int_date) %>% 
						summarise(n=n(), obs = sum(y, na.rm=TRUE), .groups='drop'), 
						by=c('round'))
		}else{
			pred = bind_cols(pred, m_d %>% group_by(round, median_int_date) %>% 
						summarise(n=n(), obs = sum(y, na.rm=TRUE), .groups='drop'))
		}
		pred = pred %>% mutate(
						label= paste(labels[i], i_mut, sep='-'),
						corr=m_output$min_qic_corr)
		write_tsv(pred, paste(out_file, '_pred.tsv', sep=''))
	}
}

