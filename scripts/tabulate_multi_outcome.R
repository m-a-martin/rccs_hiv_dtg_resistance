suppressMessages(library(tidyverse))
suppressMessages(library(argparser))
source('scripts/utils.R')


format_col = function(x1, x2){
		return(paste(paste(paste(x1, " (", sep=''),
			round((x1 / x2)*100, 1), sep=''),
			'%)', sep=''))
}

chisq_p = function(x, a, b){
	return(chisq.test(rbind(na.omit(x[[a]]), na.omit(x[[b]])))$p.value)
}


multi_outcome_table = function(hiv_dr_cat, outcomes, table_vars){
	# initiate and add overall column
	t = list()
	t[[1]] = tibble(var='overall', v=NA, outcome='n', val=nrow(hiv_dr_cat))
	for (outcome in outcomes){
		t[[length(t) + 1]] = tibble(var='overall', v=as.character(NA), 
			outcome=outcome, val=sum(hiv_dr_cat[[outcome]]))
		for (var in table_vars){
			t[[length(t) + 1]] = hiv_dr_cat[hiv_dr_cat[outcome][[1]],] %>% 
				group_by_at(var) %>% summarise(val=n()) %>%
				mutate(
					var = var,
					outcome = outcome) %>%
				rename(v := !!var) %>%
				mutate(v = as.character(v))
		}
	}
	# add overall rows for variables
	for (var in table_vars){
		t[[length(t) + 1]] = hiv_dr_cat %>% 
				group_by_at(var) %>% 
				summarise(val=n()) %>%
				mutate(
					var = var,
					outcome = 'n') %>%
				rename(v := !!var) %>%
				mutate(v = as.character(v))
	}
	# pivot
	t = bind_rows(t) %>%
		pivot_wider(names_from=outcome, values_from=val)
	# add blank rows for variables
	# todo do this above	
	for (var in table_vars){
		t = bind_rows(t, 
			tibble(var=var))
	}
	t = arrange(t, var != 'overall', var, !is.na(n))
	# add percentages 
	for (outcome in outcomes){
		t[[outcome]] = (t[c('n', outcome)] %>% mutate(k = 
			if_else(is.na(.[[1]]), NA, 
				paste(
					paste(
						format(.[[2]], big.mark=','),
						round(100*.[[2]]/.[[1]],1),
						sep=' ('),
					'%)', sep=''))))$k
	}
	return(t %>% mutate(n = format(n, big.mark=',')))
}




#### --------------- ####
#### MAIN CODE BLOCK ####
#### --------------- ####
p <- arg_parser("round prev rr table")
p <- add_argument(p, "--out", help="output file name", 
	nargs=1)
p <- add_argument(p, "--dat", help="dat file path", 
	nargs=1)
p <- add_argument(p, "--filter", help="filter string for data", 
	nargs=1, default="TRUE")
p <- add_argument(p, "--outcomes", help="data outcomes", 
	nargs=Inf, type="character")
p <- add_argument(p, "--tableVars", help="table variables", 
	nargs=Inf)
args <- parse_args(p)
print(args$outcomes)
#args$dat = 'data/rakai_drug_resistance_categorized_R15_R20.tsv'
#args$tableVars = c("round", "vl_copies_cat")
#args$filter = "!is.na(sequence_id)"
#args$outcomes = c('!is.na(nnrti)', '!is.na(nrti)')
#args$outcomeLabel = c('nnrti', 'nrti')


# read in data and filter out bad rounds and age categories
hiv_dr_cat = read_tsv(args$dat, show_col_types=FALSE)


filtered_hiv_dr_cat = hiv_dr_cat %>%
	filter(eval(parse(text=args$filter)))

# add outcome columns
for (outcome in args$outcomes){
	filtered_hiv_dr_cat = filtered_hiv_dr_cat %>%
		mutate(!!outcome := eval(parse(text=outcome)))
}
t = multi_outcome_table(filtered_hiv_dr_cat, args$outcomes, args$tableVars)

if (!dir.exists('tables')){dir.create('tables')}

write_tsv(t, 
	paste(c('tables/', args$out, '.tsv'), collapse=''), 
	na='')
