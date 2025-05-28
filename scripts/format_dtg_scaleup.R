suppressMessages(library(tidyverse))
source('scripts/utils.R')


numeric_to_datetime = function(x){
	x = as.numeric(x)
	yr = floor(x)
	yr_beg = as.POSIXct(paste0(yr, '-01-01'))
	yr_end = as.POSIXct(paste0(yr+1, '-01-01'))
	yr_length = yr_end - yr_beg
	yr_frac = x %% 1
	date = yr_beg + yr_length*yr_frac
	return(format(date, format='%Y-%m-%d'))
}


quarters = c('Q1' = '')
dtg = read_csv('data/dtg_scaleup.csv', show_col_types=FALSE) %>%
	# Did not get data from Kalisizo
	filter(clinic != "Kalisizo") %>%
	pivot_longer(-c('clinic')) %>%
	mutate(
		sex = str_split(name, "_", simplify=TRUE)[,1],
		year =  as.Date(numeric_to_datetime(as.numeric(str_split(name, "_", simplify=TRUE)[,2]) + 
			as.numeric(gsub("Q", "", str_split(name, "_", simplify=TRUE)[,3]))/4-0.125)),
		type = str_split(name, "_", simplify=TRUE)[,4],
		name = paste(
			str_split(name, "_", simplify=TRUE)[,2],
			str_split(name, "_", simplify=TRUE)[,3],
			sep = ' ')) %>%
	pivot_wider(names_from=type, values_from=value)%>%
	# drop mising data 
	filter(!is.na(DTG))

dtg = bind_rows(
	dtg, 
	dtg %>% 
		group_by(sex, name, year) %>%
		summarise(`Total ART` = sum(`Total ART`), DTG = sum(DTG), .groups='drop') %>%
		mutate(clinic = 'all')) %>%
	# add probabilities and CIs
	mutate(
		p = DTG / `Total ART`,
		s = sqrt(p*(1-p)/`Total ART`),
		lcl = p + qnorm(0.025) * s,
		lcl = if_else(lcl < 0, 0, lcl),
		ucl = p + qnorm(0.975) * s,
		ucl = if_else(ucl > 1, 1, ucl)) 

write_tsv(dtg, 'data/dtg_scaleup_format.csv')

# make a table of the overall values
dtg_table = dtg %>% filter(clinic == "all") %>%
	mutate(
		`Total ART` = format(`Total ART`, big.mark=','),
		DTG = format(DTG, big.mark=','),
		value = paste0(round(p*100,1), ' (', round(lcl*100,1), ', ', round(ucl*100,1), ')')) %>% 
	select(name, sex, `Total ART`, DTG, value) 

dtg_table = dtg_table %>% 
	filter(sex == 'men') %>%
	rename(
		`Men Total ART` = `Total ART`,
		`Men DTG` = `DTG`,
		`Men %` = value) %>%
	select(-sex) %>%
	left_join(
		dtg_table %>% filter(sex == 'women') %>%
			rename(
				`Women Total ART` = `Total ART`,
				`Women DTG` = `DTG`,
				`Women %` = value) %>%
			select(-sex),
			by=c('name'))

write_tsv(dtg_table, 'tables/dtg_by_quarter.tsv')