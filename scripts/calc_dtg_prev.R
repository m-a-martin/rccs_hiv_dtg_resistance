suppressMessages(library(tidyverse))
suppressMessages(library(geepack))
suppressMessages(library(emmeans))
suppressMessages(library(patchwork))
suppressMessages(source('scripts/utils.R'))
 

clinic_dict = c(
	"23" = "Kasasa",
	"32" = "Rakai",
	"120" = "Kamulegu",
	"25" = "Lwanda",
	"21" = "Kifamba",
	"4" = "Kasaali",
	"105" = "RHSP",
	"24" = "Mitukula",
	"121" = "Kyanamukaka",
	"12" = "Kabira",
	"89" = "Kasensero",
	"1" = "Kalisizo",
	"16" = "Kakuuto",
	"122" = "Lyantonde",
	"117" = "TASO")

# R20 PLHIV on treatment
clinic_dat = read_tsv('data/rakai_drug_resistance_categorized_R15_R20.tsv', show_col_types=FALSE) %>%
	filter(finalhiv == 'P' & round == 20 & arv == 1) %>% select('arvsourc2') %>%
	group_by(arvsourc2) %>%
	summarise(n=n(),.groups='drop') %>%
	mutate(
		clinic = clinic_dict[as.character(arvsourc2)],
		p = n / sum(n)) %>%
	arrange(-p)

nrow(clinic_dat %>% filter(!is.na(clinic)))

dtg = read_csv('data/dtg_scaleup.csv', show_col_types=FALSE) %>%
	pivot_longer(-c('clinic')) %>% 
	# Did not get data from Kalisizo
	filter(clinic != "Kalisizo") %>%
	mutate(
		sex = str_split(name, "_", simplify=TRUE)[,1],
		year =  as.numeric(str_split(name, "_", simplify=TRUE)[,2]) + 
			as.numeric(gsub("Q", "", str_split(name, "_", simplify=TRUE)[,3]))/4-0.125,
		type = str_split(name, "_", simplify=TRUE)[,4]) %>%
	select(-name) %>%
	pivot_wider(names_from=type, values_from=value) %>%
	group_by(year) %>%
	mutate(year_idx = cur_group_id()) %>%
	group_by(clinic) %>% mutate(clinic_idx = cur_group_id()) %>%
	arrange(clinic_idx, year_idx) %>%
	mutate(year = as.character(year))

# get and filter data points
nrow(dtg %>% filter(is.na(DTG)))
nrow(dtg)
dtg = dtg %>% filter(!is.na(DTG))

clinic_dat %>% inner_join(dtg %>% select(clinic) %>% unique() %>% summarise(dtg_dat = TRUE)) %>% arrange(-p)
nrow(clinic_dat %>% inner_join(dtg %>% select(clinic) %>% unique() %>% summarise(dtg_dat = TRUE)))
sum((clinic_dat %>% inner_join(dtg %>% select(clinic) %>% unique() %>% summarise(dtg_dat = TRUE)))$p)

write_tsv(
	clinic_dat %>% 
		inner_join(dtg %>% select(clinic) %>% unique() %>% summarise(dtg_dat = TRUE)) %>% arrange(clinic) %>%
		select(clinic, n, p) %>% arrange(-n) %>%
		mutate(
			p = round(100*p,1)),
	'tables/tx_by_clinic.tsv')

summary(geeglm(DTG ~ year + offset(log(`Total ART`)), 
	data=dtg %>% filter(`Total ART` > 0),
	id = clinic_idx,
	corstr='independence',
	waves=year_idx, 
	family=poisson(link="log")))

summary(geeglm(DTG ~ year + sex + offset(log(`Total ART`)), 
	data=dtg %>% filter(`Total ART` > 0),
	id = clinic_idx,
	corstr='independence',
	waves=year_idx, 
	family=poisson(link="log")))


summary(geeglm(DTG ~ year + clinic + offset(log(`Total ART`)), 
	data=dtg %>% filter(`Total ART` > 0),
	id = clinic_idx,
	corstr='independence',
	waves=year_idx, 
	family=poisson(link="log")))
