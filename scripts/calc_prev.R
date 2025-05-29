suppressMessages(library(tidyverse))
suppressMessages(library(geepack))
suppressMessages(library(emmeans))
suppressMessages(library(patchwork))
suppressMessages(source('scripts/utils.R'))
 

# todo move all of this data processing to data processing script
d = read_tsv('data/rakai_drug_resistance_categorized_R15_R20.tsv', show_col_types=FALSE)  %>%
	mutate(
		numeric_round = as.numeric(round),
		round = as.character(round))

models_to_run = read_tsv('config/models.tsv', show_col_types=FALSE, comment='#')
outcomes = models_to_run$outcome
filters = models_to_run$filter
labels = models_to_run$label
models = models_to_run$model

min_outcomes = 20
for (i in 1:length(outcomes)){
	print(i)
	print(labels[i])
	out_file = paste('models/', labels[i], sep='')
	i_d = d %>%
		# apply appropriate filter
		filter(eval(parse(text=filters[i]))) %>%
		# define outcome variable
		mutate(y = eval(parse(text=outcomes[i])))	
	# fit model by chosing the best fit correlation structure
	m_output = run_process_model(
		i_d, 
		as.formula(models[i]), 
		out_file,
		corstr=if_else(sum(i_d$y, na.rm=TRUE) > min_outcomes, 'fit', 'independence'), 
		id_col='study_id',
		order_col='round')
	# output model summary
	m_summary = as_tibble(m_output$o, rownames='var') %>% 
			mutate(
				label = labels[i],
				corr = m_output$min_qic_corr)
	# add median int dates if round in pred
	if ('round' %in% colnames(m_summary)){
		m_summary = m_summary %>% left_join(i_d %>% select(round, round_min_date, round_mid_date, round_max_date) %>% unique(), by='round')
	}
	write_tsv(m_summary, 
		paste(out_file, '.tsv', sep=''))
	# predict values in each round
	pred = as_tibble(emmeans(
				m_output$m,  
				as.formula(paste('~', gsub('y ~', '', models[i]), sep='')),
				type="response")) %>%
				rename(c(
					'fit' = 'rate',
					'se.fit' = 'SE',
					'lwr'= 'lower.CL', 
					'upr'='upper.CL'))
	# add number of observations
	pred = pred %>% left_join(
		i_d %>% 
			group_by_at(str_split(gsub(" ", "", gsub('y ~', '', models[i])), '\\+', simplify=TRUE)) %>%
			summarise(n=sum(!is.na(y)), obs = sum(y, na.rm=TRUE), .groups='drop'))
	# add median int dates if round in pred
	if ('round' %in% colnames(pred)){
		pred = pred %>% left_join(i_d %>% select(round, round_min_date, round_mid_date, round_max_date) %>% unique(), by='round')
	}
	# add labels
	pred = pred %>% mutate(
					label=labels[i],
					corr=m_output$min_qic_corr)
	write_tsv(pred, paste(out_file, '_pred.tsv', sep=''))
}


