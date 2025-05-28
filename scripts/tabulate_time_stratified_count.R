suppressMessages(library(tidyverse))
suppressMessages(library(argparser))
source('scripts/utils.R')


format_p = function(x){
	return(if_else(x < 1E-4, "<1e-4", as.character(round(x,4))))
}


time_stratified_table = function(hiv_dr_cat, table_vars){
	sum_dat = list()
	sum_dat[['oa']] = hiv_dr_cat %>% group_by(round) %>%
		summarise(val=format(n(), big.mark=',')) %>%
		mutate(s = NA, s_var = NA, round=as.character(round))
	for (i_var in table_vars){
		if (class(hiv_dr_cat[[i_var]]) != "numeric"){
			var_sum_dat = hiv_dr_cat %>%
				rename(s_var = !!i_var) %>%
				group_by(round, s_var) %>%
				summarise(n=n(), .groups='drop')
			sum_dat[[i_var]] = 
				var_sum_dat %>% group_by(round) %>%
					mutate(
						round = as.character(round),
						p = n/sum(n),
						val =  
							paste(paste(paste(format(n, big.mark=','), ' (', sep=''), round(100*p,1), sep=''), '%)', sep='')) %>%
					select(round, s_var, val) %>%
					mutate(s = i_var)
		}else{
			var_sum_dat = hiv_dr_cat %>%
				rename(s_var := !!i_var) %>%
				group_by(round) %>%
				summarise(n=paste(c(median(s_var), " [", IQR(s_var), "]"), collapse=''), .groups='drop') %>%
				mutate(s_var := !!i_var)
			sum_dat[[i_var]] = 
				var_sum_dat %>% group_by(round) %>%
					mutate(
						round = as.character(round),
						val = format(n, big.mark=',')) %>%
					select(round, s_var, val) %>%
				mutate(s = i_var)
			print(sum_dat[[i_var]])
		}
	}
	out = bind_rows(sum_dat)
	out = out %>% pivot_wider(names_from=round, values_from=val)
	return(out)
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
p <- add_argument(p, "--tableVars", help="table variables", 
	nargs=Inf)
args <- parse_args(p)


#args$dat = 'data/rakai_drug_resistance_categorized.tsv'
#args$tableVars = c("ageyrs", "age_cat", "comm_type", "sex", "vl_cat")

hiv_dr_cat = read_tsv(args$dat, show_col_types=FALSE)  %>%
	filter(eval(parse(text=args$filter)))

#hiv_dr_cat = read_tsv('data/rakai_drug_resistance_categorized_R15_R20.tsv', show_col_types=FALSE) 

t = time_stratified_table(hiv_dr_cat, args$tableVars)


if (!dir.exists('tables')){dir.create('tables')}

write_tsv(t, 
	paste(c('tables/', args$out, '.tsv'), collapse=''), 
	na='')


