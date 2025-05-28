suppressMessages(library(tidyverse))
suppressMessages(library(geepack))
suppressMessages(library(emmeans))
suppressMessages(library(patchwork))
suppressMessages(source('scripts/utils.R'))
 
dodge = c('PLHIV' = -20, 'tx' = 0, 'supp' = 20)


pred = bind_rows(
	read_tsv('models/hiv_among_par_pred.tsv', show_col_types=FALSE),
	read_tsv('models/tx_among_par_pred.tsv', show_col_types=FALSE),
	read_tsv('models/supp_among_par_pred.tsv', show_col_types=FALSE)) %>%
	mutate(label = str_split(label, '_', simplify=TRUE)[,1]) %>%
	mutate(label = ordered(label, levels=unique(label)),
		median_int_date = as.Date(median_int_date)+ dodge[label])

d = read_tsv('data/rakai_drug_resistance_categorized_R15_R20.tsv', show_col_types=FALSE)  %>%
	group_by(round) %>%
	summarise(
		min_int_date = as.Date(quantile(int_date, 0.01, na.rm=TRUE)),
		median_int_date = as.Date(median(int_date, na.rm=TRUE)),
		max_int_date = as.Date(quantile(int_date, 0.99, na.rm=TRUE)),
		.groups='drop') %>%
	merge(pred %>% select(round) %>% unique(),
		how='inner', on='round')

recs = as_tibble(d %>% slice(seq(nrow(d),1,-2)) %>%
	select(min_int_date, max_int_date) %>%
	rename(c('starts' = 'min_int_date',
		'ends' = 'max_int_date')) %>%
	mutate(starts = as.Date(starts), ends=as.Date(ends)))

ymax=30
pB = ggplot(pred) + 
	geom_rect(data=recs, aes(xmin=starts, xmax=ends, ymin=0, ymax=ymax), fill='#eaeaea') +
	geom_line(
		aes(x=median_int_date, y=fit*100, group=label),
		linewidth=1.5, color='white') + 
	geom_point(
		aes(x=median_int_date, y=fit*100, group=label),
		fill='white', color='white', shape=21, size=2) +
	geom_errorbar(
		aes(x=median_int_date, ymin=lwr*100-0.075, ymax=upr*100+0.075, group=label),
		color='white', linewidth=1.5, width=0) +
	geom_errorbar(
		aes(x=median_int_date, ymin=lwr*100, ymax=upr*100, group=label, color=label),
		width=0, linewidth=1) +
	geom_line(
		aes(x=median_int_date, y=fit*100, group=label, color=label),
		linewidth=1) + 
	geom_point(
		aes(x=median_int_date, y=fit*100, group=label, color=label, fill=label),
		size=1.5) +
	xlab('date') + 
	ylab('prevalence among\nRCCS participants (%)') +
	ggtitle('RCCS participants') +
	scale_y_continuous(limits=c(0,ymax), expand=c(0,0)) + 
	scale_color_manual(values=colors, labels=labels, name=NULL) +
	scale_fill_manual(values=colors, labels=labels, name=NULL) +
	guides(color=guide_legend(position = "inside")) +
	gtheme +
	theme(
		legend.justification.inside = c(1, 1),
		legend.position.inside=c(1,1),
		panel.border = element_blank(),
		axis.line = element_line(colour = "#333333"),
		panel.grid.major=element_blank(),
		panel.grid.minor=element_blank(),
		legend.key = element_blank()) 


#### PANEL C ####
dtg = read_tsv('data/dtg_scaleup_format.csv', show_col_types=FALSE)
recs = recs %>%
	mutate(starts = as.Date(starts), ends=as.Date(ends))

pC = ggplot(dtg %>% filter(clinic == 'all')) + 
	geom_rect(data=recs, aes(xmin=starts, xmax=ends, ymin=-2.5, ymax=102.5), fill='#eaeaea') +
	geom_line(
		aes(x=year, y=p*100, group=sex),
		linewidth=1.5, color='white') + 
	geom_point(
		aes(x=year, y=p*100, group=sex),
		fill='white', color='white', shape=21, size=2) +
	geom_errorbar(
		aes(x=year, ymin=lcl*100, ymax=ucl*100, group=sex),
		color='white', linewidth=1.5, width=0) +
	geom_errorbar(
		aes(x=year, ymin=lcl*100, ymax=ucl*100, group=sex, color=sex),
		width=0, linewidth=1) +
	geom_line(
		aes(x=year, y=p*100, group=sex, color=sex),
		linewidth=1) + 
	geom_point(
		aes(x=year, y=p*100, group=sex, color=sex),
		fill='#eaeaea', size=1.5) +
	xlab('date') + 
	ylab('DTG use (%)') +
	ggtitle('PLHIV on treatment') +
	scale_y_continuous(limits=c(-2.5, 102.5), expand=c(0,0)) + 
	scale_x_date(limits=c(min(d$median_int_date), max(d$max_int_date))) + 
	scale_color_manual(values=colors, name=NULL) +
	guides(color=guide_legend(position = "inside")) +
	gtheme +
	theme(
		legend.justification.inside = c(0, 1),
		legend.position.inside=c(0,1),
		panel.border = element_blank(),
		axis.line = element_line(colour = "#333333"),
		panel.grid.major=element_blank(),
		panel.grid.minor=element_blank(),
		legend.key = element_blank()) 


p =  (pB + pC) + 
	plot_annotation(
		tag_levels = list(c('A', 'B'))) & 
	theme(
		plot.tag = element_text(color='#333333', size=16, face='bold'))


ggsave('figures/pdf/hiv_prev_among_par.pdf', p, width=6.4*2, height=4.2)
ggsave('figures/png/hiv_prev_among_par.png', p, width=6.4*2, height=4.2)


ggsave('figures/poster_dtg_prev_among_plhiv.pdf', pC, width=4.8, height=3.6)
