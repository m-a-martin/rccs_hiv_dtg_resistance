suppressMessages(library(tidyverse))
suppressMessages(library(geepack))
suppressMessages(library(emmeans))
suppressMessages(library(patchwork))
suppressMessages(library(argparser))
suppressMessages(source('scripts/utils.R'))


#### --------------- ####
#### MAIN CODE BLOCK ####
#### --------------- ####
p <- arg_parser("risk factor table")
p <- add_argument(p, "--out", help="output file name", 
	nargs=1)
p <- add_argument(p, "--prev", help="estimated prevalence file path", 
	nargs=1)
p <- add_argument(p, "--rr", help="estimated rr file path", 
	nargs=1)
p <- add_argument(p, "--rrAdj", help="estimated rr file path", 
	nargs=1)
args <- parse_args(p)
options(scipen=999)

#args$prev = 'models/tmp_pred.tsv'
#args$rr = 'models/tmp_rr.tsv'
#args$rrAdj = 'models/tmp_rr_adj.tsv'

format_label = function(x){
	return(unlist(lapply(str_split(x, "_"), function(z){paste(z[-c(1,2,3,4)], collapse='_')})))
}


preds = read_tsv(args$prev, show_col_types=FALSE) %>% 
	select(label, var, n, obs) %>%
	mutate(
		label = format_label(label),
		p = format_digit(100*obs / n)) %>%
	unite("obs", obs:p, sep=' (') %>%
	mutate(obs = paste(obs, "%)", sep=''),
				n = as.character(n)) %>%
	pivot_longer(-c(label, var))


rr_pval = read_tsv(args$rr, show_col_types=FALSE) %>%
	filter(var != '(Intercept)') %>%
	select(label, var, RR, LCI, UCI, P) %>%
	mutate(
		P = if_else(P < 1E-4, "<0.0001", as.character(round(signif(P, 2), 4))),
		label = format_label(label)) %>%
	mutate(
		var = apply(.[c('label', 'var')], 1, function(z){gsub(z[1], "", z[2])}),
		across(RR:UCI, ~format_digit(.x, 2))) %>%
	unite("rr", RR:LCI, sep=' (') %>%
	unite("rr", rr:UCI, sep=', ') %>%
	mutate(rr = paste(rr, ")", sep='')) %>%
	pivot_longer(-c(label, var))


adj_rr_pval = read_tsv(args$rrAdj, show_col_types=FALSE) %>%
	filter(!grepl('round', var), var != '(Intercept)') %>%
	select(label, var, RR, LCI, UCI, P) %>%
	mutate(
		P = if_else(P < 1E-4, "<0.0001", as.character(round(signif(P, 2), 4))),
		label = gsub("_round", "", format_label(label))) %>%
	mutate(
		var = apply(.[,c('label', 'var')], 1, function(z){gsub(z[1], "", z[2])}),
		across(RR:UCI, ~format_digit(.x, 2))) %>%
	unite("rr", RR:LCI, sep=' (') %>%
	unite("rr", rr:UCI, sep=', ') %>%
	mutate(rr = paste(rr, ")", sep='')) %>%
	rename(rr_adj = rr, P_adj = P) %>%
	pivot_longer(-c(label, var))


table = bind_rows(preds, rr_pval, adj_rr_pval) %>%
	pivot_wider(names_from=name, values_from=value, values_fill='ref') 

table = bind_rows(table, table %>% select(label) %>% unique() %>% mutate(var=NA))

table = table %>%
	arrange(label, !is.na(var), var)

write_tsv(table, paste(c('tables/', args$out, ".tsv"), collapse=''), na="")


