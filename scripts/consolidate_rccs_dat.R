suppressMessages(require(tidyverse))
suppressMessages(require(haven))
suppressMessages(require(readxl))


#### 1. DEFINE RCCS VARIABLES TO GET  ####
rccs_vars = c('study_id', 'round', 'int_date', 'sex', 'ageyrs', 'age_cat', 'comm_cat', 'comm_type', 
		'ever_arv', 'arv',
		'arvsourc2',
		#'last_non_arv', 'last_non_arv_date', 
		'first_arv', 'first_arv_date', 
		'finalhiv',
		'first_pos', 'first_pos_date',
		'last_neg', 'last_neg_date', 
		'first_supr', 'first_supr_date',
		#'last_non_supr', 'last_non_supr_date',
		'pre_treatment',
		#'est_tx_date', 'est_tx_duration',
		'copies', 'numeric_copies', 'vl_copies_cat')

arv_vars = c('cuarvmed', 'arvmed', 'ltemmed', 
			'hivcare', 'hivrslt', 'rhivrlst', 'hivquest', 'rhivever')

#### 2. READ IN RAKAI DATA  ####
rccs_data_file = 'data/R20_15.csv'
#rccs_backup = read_csv(rccs_data_file) %>%
#	mutate(
#		round = as.numeric(str_replace(str_remove(round, "^R0+"), 'S', '.1')))
load("data/R20_15.RData")
#save(rccs_backup, file="data/R20_15.RData")
# fitler for eligibility
# age 15-49
# locate = 1,2,3,4
# have hiv status
# not recruited as part of HTR
# exclude out migrants NEED TO ASK KATE WHICH VARIABLE TO USE 
rccs = rccs_backup %>% 
	filter(
		ageyrs >= 15 & ageyrs <= 49 &
		locate %in% c(1,2,3,4) & 
		finalhiv %in% c("N", "P") & 
		(htr_aim == 8 | is.na(htr_aim) | (htr_aim == 2 & comm_num != 600)) &
		(vicinity ==8 | is.na(vicinity))) %>%
		select(any_of(c(rccs_vars, arv_vars))) %>%
	mutate(sex = str_to_upper(sex))

# fill in missing R15 community types
rccs = rccs %>% 
	left_join(
		read_stata('data/RCCSdata_R001_R019_VOIs.dta') %>% select(study_id, round, comm_type) %>%
		mutate(
			round = as.numeric(str_replace(str_remove(round, "^R0+"), 'S', '.1'))) %>% 
		mutate(
			comm_type2 = setNames(names(attributes(comm_type)$labels), 
				as.character(attributes(comm_type)$labels))[as.character(comm_type)]) %>%
			mutate(comm_type2 = case_when(
				comm_type2 == "Trading" ~ "Trading community",
				comm_type2 == "Agrarian" ~ "Agrarian community",
				comm_type2 == "Fishing" ~ "Fish landing site")) %>%
			select(-comm_type),
		by=c('study_id', 'round')) %>%
	mutate(
		comm_type = case_when(
			is.na(comm_type) ~ comm_type2,
			!is.na(comm_type) ~ comm_type)) %>%
	select(-comm_type2)

# read in dates
rccs_dates_file = 'data/r20_15_dates.xlsx'
rccs_dates = read_excel(rccs_dates_file) %>%
	mutate(
		round = as.numeric(str_replace(str_remove(round, "^R0+"), 'S', '.1'))) %>%
	filter(!is.na(int_date)) %>%
	select(-curr_id) 

# merge
rccs_merged = rccs_dates %>% right_join(rccs, by=c('study_id', 'round'))
stopifnot(max((rccs_merged %>% group_by(study_id, round) %>% summarise(n=n(), .groups='drop'))$n) == 1)


#### 4. ADD NUMERIC_COPIES ####
rccs_merged  = rccs_merged %>%
	mutate(
		numeric_copies = as.numeric(case_when(
			finalhiv == 'N' ~ NA,
			copies == 'QNS' ~ NA,
			copies == 'INV.IC' ~ NA,
			copies == "" ~ NA,
	 		grepl("<", copies) ~ "-1",
			grepl("Not Det", copies) ~ "-1",
			as.numeric(copies) < 40 ~ "-1",
			copies == "BD" ~ "-1",
			copies == "ND" ~ "-1",
			TRUE ~ gsub(",", "", copies))),
		log10_copies = case_when(
			numeric_copies == -1 ~ -1,
			numeric_copies == 0 ~ 0,
			TRUE ~ log10(numeric_copies)),
		vl_copies_cat = cut(log10_copies, c(-Inf, -1, 0, 3, 4, 5, Inf)))


#### 3. GENERATE COMPOSITE TREATMENT VARIABLES ####
treatment_data =
	bind_rows(
		rccs_merged %>% 
			select(study_id, round, int_date, finalhiv, cuarvmed, arvmed, ltemmed, 
				hivcare, hivrslt, rhivrlst, hivquest, rhivever),
		read_stata('data/RCCSdata_R001_R019_VOIs.dta') %>%
			mutate(round = as.numeric(str_replace(str_remove(round, "^R0+"), 'S', '.1'))) %>%
			filter(round < 15) %>%
			select(any_of(c('study_id', 'round', 'cuarvmed', 'arvmed', 'ltemmed', 
				'hivcare', 'hivrslt', 'rhivrlst', 'hivquest', 'rhivever'))) %>%
			left_join(read_excel('data/RCCSdata_rhiv.xlsx') %>%
					mutate(round = as.numeric(str_replace(str_remove(round, "^R0+"), 'S', '.1'))),
				by=c('study_id', 'round')))

# below intentionally verbose
arv = treatment_data %>% mutate(arv = case_when(
	# if not HIV positive then not on ARV
	finalhiv != "P" ~ FALSE,
	# here treating cuarvmed and arvmed both as current 
	# cuarvmed = currently on ART
		# 1 = "yes"
		cuarvmed == 1 ~ TRUE,
		cuarvmed == 2 ~ FALSE,
	# arvmed = have you ever been on ART
		# 1 = "yes"
		arvmed == 1 ~ TRUE,
		# 2 = "no"
		arvmed == 2 ~ FALSE,
	# ltemmed = are you on any long term medications
		# 2 = "no"
		# if no long term medication then not on ART
		ltemmed == 2 ~ FALSE,
	# hivcare = have you ever been to a clinic to receive HIV care
		# 2 = "no"
		# if never received HIV care then not on ART
		hivcare == 2 ~ FALSE,
	# hivrlst = what was the last test result
		# 1 = negative
		# if last test result negative then not on ART
		hivrslt == 1 ~ FALSE,
		# if last test result indeterminant then not on ART
		# 3 = indeterminant
		hivrslt == 3 ~ FALSE,
	# rhivrlst = did you receive the results
		# 2 == "no"
		# if did not receive results then not on ART
		rhivrlst == 2 ~ FALSE,
	# hivquest = have you ever requested your HIV results (R15 only)
		# 2 = "gave blood but did not request results"
		hivquest == 2 ~ FALSE,
		# 8 = "never tested"
		hivquest == 8 ~ FALSE,
	# rhivever = have you ever received hiv results
		# 2 = "no"
		rhivever == 2 ~ FALSE,
		# 8 = "never tested",
		rhivever == 8 ~ FALSE))

# do a cursory check for incompatible patterns
# taking arvs currently but report never having taken arts
rccs %>% filter(arvmed != 1 & cuarvmed == 1)
# never tested or never received results yet report having received a result
rccs %>% filter((hivquest == 2 | hivquest == 8) & rhivrlst == 1)
# never tested or never received results yet report having received a positive result
rccs %>% filter((hivquest == 2 | hivquest == 8) & hivrslt != 8)
# never received results yet report having received a result
rccs %>% filter((rhivever == 2 | rhivever == 8) & rhivrlst == 1)
# never received results yet report having received a positive result
rccs %>% filter((rhivever == 2 | rhivever == 8) & hivrslt != 8)
# rhivever & hivquest not the same
rccs %>% filter(!is.na(rhivever) & !is.na(hivquest))

# count missing data
print(paste(nrow(arv %>% filter(round >= 15) %>% filter(is.na(arv))),
	' participant-visits in R15-20 have missing data on ART status. Assuming pre_treatment.', sep=''))

# merge arv with all data
first_arv = arv %>% 
	filter(finalhiv == "P" & arv == TRUE) %>% 
	group_by(study_id) %>% 
	filter(round == min(round)) %>%
	mutate(
		first_arv = round, 
		first_arv_date = int_date
		) %>%
	select(
		study_id, 
		first_arv, 
		first_arv_date
		)

last_non_arv = arv %>%
	left_join(first_arv, by=c('study_id')) %>%
	# only get rounds before first arv round
	filter(finalhiv == "P" & arv == FALSE & round < first_arv) %>%
	group_by(study_id) %>%
	filter(round == max(round)) %>%
	mutate(last_non_arv = round, last_non_arv_date = int_date) %>%
					select(study_id, last_non_arv, last_non_arv_date)

rccs_merged = rccs_merged %>%
	left_join(
		arv %>% 
			filter(finalhiv == "P") %>%
			select(study_id, round, arv) %>% 
			left_join(first_arv,
				by=c('study_id')) %>%
			left_join(
				last_non_arv,
				by=c('study_id')) %>%
			mutate(
				first_arv = replace_na(first_arv, Inf))) %>%
	mutate(arv = if_else(finalhiv == "N", FALSE, arv))


#### 4. GENERATE FIRST_HIV STATUS VARIABLE ####
hiv =
	bind_rows(
		rccs_merged %>% 
			select(study_id, round, int_date, finalhiv),
		read_stata('data/RCCSdata_R001_R019_VOIs.dta') %>%
			mutate(
				round = as.numeric(str_replace(str_remove(round, "^R0+"), 'S', '.1')),
				int_date = as_date(int_date, format='%d/%m/%Y')) %>%
			filter(round < 15) %>%
			select(any_of(c('study_id', 'round', 'int_date', 'finalhiv'))))

first_last_hiv = hiv %>% 
	filter(finalhiv == 'N') %>%
	group_by(study_id) %>%
	filter(round == max(round)) %>%
	mutate(last_neg = round, last_neg_date = int_date) %>%
	select(study_id, last_neg, last_neg_date) %>%
	full_join(
		hiv %>%
			filter(finalhiv == 'P') %>%
			group_by(study_id) %>%
			filter(round == min(round)) %>%
			mutate(first_pos = round, 
				first_pos_date=int_date
				) %>%
			select(study_id, first_pos, 
				first_pos_date
				),
		by=c('study_id'))

rccs_merged = rccs_merged %>%
	left_join(first_last_hiv, 
		by=c('study_id'))
	

#### 5. GENERATE SUPPRESSION STATUS VARIABLE ####
# PARSING FAILURES ARE IN HERE
# TODO INVESTIGATE
supr = bind_rows(
		rccs_merged %>% 
			select(study_id, round, int_date, finalhiv, numeric_copies),
		read_stata('data/RCCSdata_R001_R019_VOIs.dta') %>%
			mutate(
				round = as.numeric(str_replace(str_remove(round, "^R0+"), 'S', '.1')),
				int_date = as_date(int_date, format='%d/%m/%Y')) %>%
			filter(round < 15) %>%
			select(any_of(c('study_id', 'round', 'int_date', 'finalhiv', 'copies')))) %>%
	filter(finalhiv == "P" & !is.na(numeric_copies))

# numeric_copies = -1 is code for BD or similar (e.g. "<150")
first_supr = supr %>% 
	filter(numeric_copies < 1000) %>%
	group_by(study_id) %>%
	filter(round == min(as.numeric(round))) %>%
	mutate(
		first_supr = as.numeric(round),
		first_supr_date = int_date) %>%
	select(study_id, first_supr, 
		first_supr_date)

# last non suppressed before first suppressed
last_non_supr = supr %>%
	left_join(first_supr, by='study_id') %>%
	mutate(first_supr = replace_na(first_supr, Inf)) %>%
	filter(round < first_supr & numeric_copies >= 1000) %>%
	group_by(study_id) %>%
	filter(round == min(as.numeric(round))) %>%
	mutate(last_non_supr = as.numeric(round), last_non_supr_date = int_date) %>%
	select(study_id, last_non_supr, last_non_supr_date)

rccs_merged = rccs_merged %>%
	left_join(first_supr, 
		by=c('study_id')) %>%
	mutate(first_supr = replace_na(first_supr, Inf)) %>%
	left_join(last_non_supr,
		by=c('study_id'))
	

#### 6. ADD PRE-TREATMENT VARIABLE ####
 rccs_merged = rccs_merged %>% 
	mutate(
		pre_treatment = case_when(
		finalhiv == "N" ~ NA,
		finalhiv == "P" & (round < first_arv & round < first_supr) ~ TRUE,
		finalhiv == "P" ~ FALSE)) 

# okay so if 
# 
#rccs_merged = rccs_merged %>%
#	mutate(
#		est_tx_date = 
#			case_when(
#				# if not HIV then NA 
#				finalhiv != "P" ~ NA,
#				# if pre_treatment then NA
#				finalhiv == "P" & pre_treatment == TRUE ~ 
#					NA,
#				# if we don't have a last non arv 
#				finalhiv == "P" & pre_treatment == FALSE & is.na(last_non_arv) ~ 
#					NA,
#				# if first suppression is before last non arv then we don't have an accurate last non arv
#				finalhiv == "P" & pre_treatment == FALSE & first_supr < last_non_arv ~ 
#					NA,
#				# suppressed viral load before first arv, so take date of first suppression
#				finalhiv == "P" & pre_treatment == FALSE & first_supr < first_arv ~ 
#					as.Date(last_non_arv_date) + (first_supr_date - last_non_arv_date)/2,
#				# report being on treatment so take date of first reported treatment
#				finalhiv == "P" & pre_treatment == FALSE & first_arv < Inf ~ 
#					as.Date(last_non_arv_date) + (first_arv_date - last_non_arv_date)/2),
#		est_tx_duration = as.numeric(as.Date(int_date) - est_tx_date))

#### 7. ADD AGE CATEGORY VARIABLE ####
breaks = c(14,24,34,49)
rccs_merged = rccs_merged %>%
	mutate(age_cat = cut(ageyrs, breaks=breaks))

#### 8. ADD COMM CATEGORY VARIABLE ####
rccs_merged = rccs_merged %>%
	mutate(comm_cat = if_else(comm_type == "Fish landing site", 'fishing', 'inland'))

write_csv(rccs_merged %>% select(any_of(c(rccs_vars, 'min_int_date', 'median_int_date', 'max_int_date'))), 
	gsub('.csv', '_consolidated.csv', rccs_data_file))
