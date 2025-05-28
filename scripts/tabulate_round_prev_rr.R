suppressMessages(library(tidyverse))
suppressMessages(library(argparser))
source('scripts/utils.R')


label_round_prev_rr_table = function(pred_file, rr_file, label_order){
	sort_cols = function(x, label_order){
		type_order = c("n" = 1, "obs" = 2, "Prev. % (95% CI)" = 3, "Prev. ratio (95% CI)" = 4, "p-value" = 5, "spacer" = 6)
		non_round_x = x[x!= "round"]
		non_round_x_split = str_split(non_round_x, "_")
		s1 = label_order[unlist(lapply(non_round_x_split, function(x){paste(x[1:(length(x)-1)], collapse='_')}))]
		s2 = type_order[unlist(lapply(non_round_x_split, function(x){paste(x[(length(x))], collapse='_')}))]
		scores = tibble(val = non_round_x, s1 = s1, s2 = s2) %>%
			arrange(s1, s2)
		return(c("round", scores$val))
	}
	all_pred = read_tsv(pred_file, show_col_types=FALSE)
	# format for table
	pred = all_pred %>%
		select(-any_of(c("n_outcome_obs", "n", "obs"))) %>%
		select(-se.fit, -df) %>%
		mutate(across(fit:upr, ~format_digit(.x*100))) %>%
		unite("val", fit:lwr, sep=' (') %>%
		unite("val", val:upr, sep=', ') %>%
		mutate(
			val = paste(val, ')', sep=''),
			type = paste(label, "Prev. % (95% CI)", sep='_')) %>%
		select(-any_of(c('label',  'corr', 'class')))
	if (any(colnames(all_pred) == 'n')){
		pred = bind_rows(
			pred,
			all_pred %>% 
				mutate(type = paste(label, "n", sep='_'), val=format(n, big.mark=',')) %>% 
				select(round, val,type),
			all_pred %>% 
				mutate(type = paste(label, "obs", sep='_'), val=format_col(obs,n)) %>% 
				select(round, val,type))
	}
	rr = read_tsv(rr_file, show_col_types=FALSE) %>%
		select(-any_of(c("value", "n_outcome_obs"))) %>%
		mutate(across(RR:UCI, ~format_digit(.x, 2))) %>%
		unite("val", RR:LCI, sep=' (') %>%
		unite("val", val: UCI, sep=', ') %>%
		group_by(label) %>%
		mutate(
			val = paste(val, ')', sep=''),
			val = if_else(`var` == "(Intercept)", 'ref', val),
			P = if_else(P < 1E-4, "<0.0001", as.character(signif(P, 2))),
			P = if_else(`var` == "(Intercept)", 'ref', P),
			`var` = as.numeric(str_split(`var`, 'round', simplify=T)[,2]),
			`var` = if_else(is.na(`var`), 
				min(`var`, na.rm=TRUE)-1,
				`var`),
			type = paste(label, "Prev. ratio (95% CI)", sep='_')) %>%
		ungroup() %>%
		rename(c("round"='var')) %>%
		select(-any_of(c('label',  'corr')))
	rr = 
		bind_rows(
			rr %>% select(-P),
			rr %>% select(-val) %>% rename(c("val" = "P")) %>%
				mutate(type = gsub('Prev.\ ratio \\(95%\ CI\\)', 'p-value', type)))
	t = bind_rows(pred, rr)
	t = t %>% select(-any_of(c('min_int_date', 'max_int_date', 'median_int_date', 'label', 'class'))) %>%
		pivot_wider(names_from=type, values_from=val)
	# add blank columns
	for (label in unique(read_tsv(pred_file, show_col_types=FALSE)$label)){
		t[paste(label, 'spacer', sep='_')] = " "
	}
	# sort
	t = t %>% select(sort_cols(colnames(t), label_order))
	if ('s' %in% colnames(t)){
		t = t %>% arrange(s, !is.na(s_var), s_var, round)
	}
	# %>%
	#	arrange(s, round, !is.na(s_var), s_var)
	return(t)
}


#### --------------- ####
#### MAIN CODE BLOCK ####
#### --------------- ####
p <- arg_parser("round prev rr table")
p <- add_argument(p, "--out", help="output file name", 
	nargs=1)
p <- add_argument(p, "--prev", help="estimated prevalence file path", 
	nargs=1)
p <- add_argument(p, "--rr", help="estimated rr file path", 
	nargs=1)
p <- add_argument(p, "--labelOrder", help="estimated rr file path", 
	nargs=Inf)
p <- add_argument(p, "--nOut", help="estimated rr file path", 
	nargs=1, default="all")
p <- add_argument(p, "--endRound", help="used to order observations to sort by order", 
	nargs=1, type="numeric")
args <- parse_args(p)
options(scipen=999)
#args$prev = 'models/tmp_prev_pred.tsv'
#args$rr = 'models/tmp_rr.tsv'
#args$labelOrder = c('viremic_tx_inT97A', 'viremic_tx_rtK103N')
#args$rr = 'models/dr_amongPar_rr.tsv'
#args$prev = 'models/dr_amongPar_prev_pred.tsv'
#args$label_order =c('plhiv', 'v_plhiv', 'pt_v_plhiv', 't_v_plhiv')
#args$labelOrder = c('nnrti', 'nrti', 'pi', 'nnrti_nrti', 'nnrti_pi', 'nrti_pi', 'nnrti_nrti_pi')
#args$labelOrder = c('nnrti', 'nrti', 'pi')

if (all(is.na(args$labelOrder))){
	prev = read_tsv(args$prev, show_col_types=FALSE) %>% filter(round == args$endRound) %>%
		arrange(-fit)
	label_order = setNames(seq(1,nrow(prev)), prev$label)
}else{
	label_order = setNames(seq(1,length(args$labelOrder)), args$labelOrder)
}


t = label_round_prev_rr_table(args$prev, args$rr, label_order)

if (args$nOut == 'one'){
	t = t %>% 
		select(-any_of(paste(names(label_order)[2:length(label_order)], '_n', sep='')))
}
if (!dir.exists('tables')){dir.create('tables')}

write_tsv(t, 
	paste(c('tables/', args$out, '.tsv'), collapse=''), 
	na='')
