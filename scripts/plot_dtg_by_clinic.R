suppressMessages(library(tidyverse))
suppressMessages(library(patchwork))
suppressMessages(source('scripts/utils.R'))
 
dodge = c('PLHIV' = -20, 'tx' = 0, 'supp' = 20)


dtg = read_csv('data/dtg_scaleup.csv', show_col_types=FALSE) %>%
	pivot_longer(-c('clinic')) %>%
	# Did not get data from Kalisizo
	filter(clinic != "Kalisizo") %>%
	mutate(
		sex = str_split(name, "_", simplify=TRUE)[,1],
		year =  as.numeric(str_split(name, "_", simplify=TRUE)[,2]) + 
			as.numeric(gsub("Q", "", str_split(name, "_", simplify=TRUE)[,3]))/4-0.125,
		type = str_split(name, "_", simplify=TRUE)[,4]) %>%
	filter(year <= 2023.125) %>%
	select(-name) %>%
	pivot_wider(names_from=type, values_from=value) %>%
	group_by(clinic, year) %>%
	summarise(`Total ART` = sum(`Total ART`), DTG = sum(DTG), p = DTG/`Total ART`) %>%
	mutate(clinic = if_else(clinic == "TASO", "TASO Masaka", clinic)) %>%
	mutate(clinic_type = ordered(
		case_when(
		clinic == "TASO Masaka" | clinic == "Lwanda" | clinic == "Kifamba" ~ clinic,
		TRUE ~ "Other"),
		levels=c("Kifamba", "Lwanda", "TASO Masaka", "Other"))) %>%
	filter(!is.na(DTG))


p = ggplot() + 
	geom_line(data=dtg,
			aes(x=year, y=p*100, group=clinic),
			linewidth=2, color='white') + 
	geom_line(data=dtg %>% filter(clinic_type == "Other"),
		aes(x=year, y=p*100, group=clinic, color=clinic_type),
		linewidth=1.5) +
	geom_line(data=dtg %>% filter(clinic_type != "Other"),
		aes(x=year, y=p*100, group=clinic, color=clinic_type),
		linewidth=1.5) +
	scale_color_manual(values=c("TASO Masaka"="#4d8322", "Lwanda"="#224d83", "Kifamba" ="#83224d", "Other" = "#848484"),
	breaks=levels(dtg$clinic_type), name=NULL) +
	xlab('date') +
	ylab('DTG use (%)') +
	guides(color=guide_legend(position = "inside")) +
	gtheme + 
	theme(
		legend.justification.inside = c(0, 1),
		legend.position.inside=c(0,1),
		panel.grid.major.x = element_blank(),
	panel.grid.minor.x = element_blank(),
panel.border = element_blank(),
		axis.line = element_line(colour = "#333333"))


ggsave('figures/pdf/dtg_clinic_among_par.pdf', p, width=6.6, height=4.4)
ggsave('figures/png/dtg_clinic_among_par.png', p, width=6.6, height=4.4)


