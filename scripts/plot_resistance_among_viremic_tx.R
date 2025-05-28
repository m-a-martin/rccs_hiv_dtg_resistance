suppressMessages(library(tidyverse))
suppressMessages(library(geepack))
suppressMessages(library(emmeans))
suppressMessages(library(patchwork))
suppressMessages(source('scripts/utils.R'))
 

## PANEL A ##
dodge = c('any'=-40, 'insti' = -20, 'nnrti'=0, "nrti" = 20, "pi" = 40)
pred = bind_rows(
	read_tsv('models/any_among_viremic_tx_pred.tsv', show_col_types=FALSE),
	read_tsv('models/insti_among_viremic_tx_pred.tsv', show_col_types=FALSE),
	read_tsv('models/nnrti_among_viremic_tx_pred.tsv', show_col_types=FALSE),
	read_tsv('models/nrti_among_viremic_tx_pred.tsv', show_col_types=FALSE),
	read_tsv('models/pi_among_viremic_tx_pred.tsv', show_col_types=FALSE)) %>%
	mutate(label = str_split(label, '_', simplify=TRUE)[,1]) %>%
	mutate(label = ordered(label, levels=unique(label)),
		median_int_date = as.Date(median_int_date)+ dodge[label]) %>%
	mutate(lwr = if_else(obs == 0, NA, lwr),
		upr = if_else(obs == 0, NA, upr))

# generate shading rectangles
# todo move this to calc_prev file
d = read_tsv('data/rakai_drug_resistance_categorized_R15_R20.tsv', show_col_types=FALSE)  %>%
	group_by(round) %>%
	summarise(
		min_int_date = as.Date(quantile(int_date, 0.01, na.rm=TRUE)),
		median_int_date = as.Date(median(int_date, na.rm=TRUE)),
		max_int_date = as.Date(quantile(int_date, 0.99, na.rm=TRUE)),
		.groups='drop') %>%
	merge(pred %>% select(round) %>% unique(),
		how='inner', on='round')

recs = as_tibble(d %>% slice(seq(5,1,-2)) %>%
	select(min_int_date, max_int_date) %>%
	rename(c('starts' = 'min_int_date',
		'ends' = 'max_int_date')) %>%
	mutate(starts = as.Date(starts), ends=as.Date(ends)))

ymax=5*ceiling(max(pred$upr, na.rm=TRUE)*100/5)
pA = ggplot(pred) + 
	geom_rect(data=recs, aes(xmin=starts, xmax=ends, ymin=-3, ymax=ymax), fill='#eaeaea') +
	geom_line(
		aes(x=median_int_date, y=fit*100, group=label),
		linewidth=2, color='white') + 
	geom_point(
		aes(x=median_int_date, y=fit*100, group=label, shape=label),
		fill='white', color='white', size=5) +
	geom_errorbar(
		aes(x=median_int_date, ymin=lwr*100, ymax=upr*100, group=label),
		color='white', linewidth=2, width=0) +
	geom_errorbar(
		aes(x=median_int_date, ymin=lwr*100, ymax=upr*100, group=label, color=label),
		width=0, linewidth=1.5) +
	geom_line(
		aes(x=median_int_date, y=fit*100, group=label, color=label),
		linewidth=1.5) + 
	geom_point(
		aes(x=median_int_date, y=fit*100, group=label, color=label, shape=label),
		fill='#eaeaea', size=3, stroke=1.5) +
	scale_y_continuous(expand=c(0,0), limits=c(-3, ymax)) +
	ylab('prevalence (%)') +
	xlab('date') + 
	scale_shape_manual(values=shapes, name=NULL) +
	scale_color_manual(values=colors, name=NULL) +
	guides(
		shape=guide_legend(position="inside"),
		color=guide_legend(position="inside")) +
	gtheme +
	theme(
		legend.justification.inside = c(1, 1),
		legend.position.inside=c(1,1),
		panel.border = element_blank(),
		axis.line = element_line(colour = "#333333"),
		panel.grid.major=element_blank(),
		panel.grid.minor=element_blank(),
		legend.key = element_blank()) 

## PANEL B ##
# any resistance risk differences
m = readRDS('models/any_among_viremic_tx.rds')
m_pairs = as_tibble(pairs(regrid(emmeans(m, "round"), "response"), infer=c(TRUE, TRUE))) %>%
	filter(
		contrast == 'round16 - round17' |
		contrast == 'round17 - round18' | 
		contrast == 'round18 - round19' | 
		contrast == 'round19 - round20') %>%
	select(contrast, estimate, lower.CL, upper.CL) %>%
	mutate(
		estimate = -1*estimate, lwr = -1*upper.CL, upr = -1*lower.CL,
		label = case_when(
			contrast == 'round16 - round17' ~ '2015 v. 2014',
			contrast == 'round17 - round18' ~ '2017 v. 2015',
			contrast == 'round18 - round19' ~ '2019 v. 2017',
			contrast == 'round19 - round20' ~ '2022 v. 2019'))

# need to get number of years separating rounds
years = read_tsv('models/any_among_viremic_tx_pred.tsv', show_col_types=FALSE) %>%
	select(round, median_int_date)
m_pairs = m_pairs %>%
	mutate(
		round1 = as.numeric(gsub("round", "", str_split(contrast, " - ", simplify=TRUE)[,1])),
		round2 = as.numeric(gsub("round", "", str_split(contrast, " - ", simplify=TRUE)[,2]))) %>%
	left_join(years %>% rename(c('round1' = 'round', 'median_int_date1' = 'median_int_date')), by='round1') %>%
	left_join(years %>% rename(c('round2' = 'round', 'median_int_date2' = 'median_int_date')), by='round2') %>%
	mutate(
		dt = as.numeric((median_int_date2 - median_int_date1)/365.25),
		estimate = estimate/dt, 
		lwr = lwr/dt,
		upr = upr/dt)

pB = ggplot(m_pairs, aes(x=label, y=estimate*100)) + 
	geom_hline(aes(yintercept=0.0), linetype='dashed', color='#4d4d4d') +
	geom_point(
		fill='white', color='white', size=5) +
	geom_errorbar(aes(ymin=lwr*100, ymax=upr*100),
		color='white', linewidth=2, width=0) +
	geom_errorbar(aes(ymin=lwr*100, ymax=upr*100), width=0, linewidth=1.5, color='#333333') +
	geom_point(color='#333333', fill='#eaeaea', shape=21, size=3, stroke=1.5) +
	ylab( expression(paste(Delta, " prevalence (%) / year"))) +
	xlab(NULL) +
	ylim(-12.5, 12.5) +
	gtheme +
	theme(panel.grid.major.x=element_blank(),
		panel.grid.minor.x=element_blank(),
		axis.text.x = element_text(size=10),
		panel.border = element_blank(),
		axis.line = element_line(colour = "#333333"))


## PANEL C ##
dodge = c('single-class'=-10, 'multi-class'=10)

pred = bind_rows(
	read_tsv('models/multi-class_among_viremic_tx_pred.tsv', show_col_types=FALSE),
	read_tsv('models/single-class_among_viremic_tx_pred.tsv', show_col_types=FALSE)) %>%
	mutate(label = str_split(label, '_', simplify=TRUE)[,1]) %>%
	mutate(label = ordered(label, levels=unique(label)),
		median_int_date = as.Date(median_int_date)+ dodge[label]) %>%
	mutate(lwr = if_else(obs == 0, NA, lwr),
		upr = if_else(obs == 0, NA, upr))


ymax=5*ceiling(max(pred$upr, na.rm=TRUE)*100/5)
pC = ggplot(pred) + 
	geom_rect(data=recs, aes(xmin=starts, xmax=ends, ymin=-3, ymax=ymax), fill='#eaeaea') +
	geom_line(
		aes(x=median_int_date, y=fit*100, group=label),
		linewidth=2, color='white') + 
	geom_point(
		aes(x=median_int_date, y=fit*100, group=label),
		fill='white', color='white', size=5, shape=21) +
	geom_errorbar(
		aes(x=median_int_date, ymin=lwr*100, ymax=upr*100, group=label),
		color='white', linewidth=2, width=0) +
	geom_errorbar(
		aes(x=median_int_date, ymin=lwr*100, ymax=upr*100, group=label, color=label),
		width=0, linewidth=1.5) +
	geom_line(
		aes(x=median_int_date, y=fit*100, group=label, color=label),
		linewidth=1.5) + 
	geom_point(
		aes(x=median_int_date, y=fit*100, group=label, color=label),
		fill='#eaeaea', size=3, stroke=1.5, shape=21) +
	scale_y_continuous(expand=c(0,0), limits=c(-3, ymax)) +
	ylab('prevalence (%)') +
	xlab('date') + 
	scale_color_manual(values=colors, name=NULL) +
	guides(color=guide_legend(position="inside")) +
	gtheme +
	theme(
		legend.justification.inside = c(1, 1),
		legend.position.inside=c(1,1),
		panel.border = element_blank(),
		axis.line = element_line(colour = "#333333"),
		panel.grid.major=element_blank(),
		panel.grid.minor=element_blank(),
		legend.key = element_blank()) 


#### PANEL D ####
files = Sys.glob(file.path('models/mutations', "*tx*_pred.tsv"))
mut_pred = do.call('rbind',lapply(files,read_tsv)) %>%
	filter(round == 20) %>%
	arrange(-fit) %>% 
	mutate(mut = str_split(label, "-", simplify=TRUE)[,3]) %>%
	left_join(
		read_tsv('config/mut_class.tsv', show_col_types=FALSE),
		by='mut') %>%
	mutate(mut = ordered(mut, levels=mut))

pD = ggplot(mut_pred) + 
	geom_point(
		aes(x=mut, y=fit*100, group=label),
		fill='white', color='white', size=5, shape=21) +
	geom_errorbar(
		aes(x=mut, ymin=lwr*100, ymax=upr*100, group=label),
		color='white', linewidth=2, width=0) +
	geom_errorbar(
		aes(x=mut, ymin=lwr*100, ymax=upr*100, group=label, color=class),
		width=0, linewidth=1.5) +
	geom_point(
		aes(x=mut, y=fit*100, group=label, color=class),
		fill='#eaeaea', size=3, stroke=1.5, shape=21) +
	ylab('prevalence (%)') +
	xlab(NULL) + 
	ylim(0,5*ceiling(max(mut_pred$upr)*100/5)) +
	scale_color_manual(values=colors, name=NULL) +
	guides(color=guide_legend(position="inside")) +
	gtheme +
	theme(
		legend.justification.inside = c(1, 1),
		legend.position.inside=c(1,1),
		panel.border = element_blank(),
		axis.line = element_line(colour = "#333333"),
		panel.grid.major.x=element_blank(),
		panel.grid.minor.x=element_blank(),
		legend.key = element_blank(),
		axis.text.x = element_text(angle=90))

p = (pA + pB + plot_layout(widths = c(1.25, 1))) / (pC + pD + plot_layout(widths = c(1.25, 1))) + 
	plot_annotation(
		tag_levels = list(c('A', 'B', 'C', 'D')),
		title = 'treatment-experienced viremic PLHIV') &
	theme(
		plot.tag = element_text(color='#333333', size=16, face='bold'),
		plot.title = element_text(color='#333333', size=22, hjust = 0.5))

ggsave('figures/pdf/resistance_among_viremic_tx.pdf', p, width=6.4*1.8, height=4.8*1.8)
ggsave('figures/png/resistance_among_viremic_tx.png', p, width=6.4*1.8, height=4.8*1.8)


p = (pA + xlab(NULL)) + plot_spacer()  + pD + plot_layout(widths = c(1.5, 0.25, 1))+ 
	plot_annotation(
		title = 'treatment-experienced viremic PLHIV')&
	theme(
		plot.tag = element_text(color='#333333', size=16, face='bold'),
		plot.title = element_text(color='#333333', size=22, hjust = 0.5))
ggsave('figures/poster_resistance_among_viremic_tx.pdf', p, width=10.8, height=4.2)

