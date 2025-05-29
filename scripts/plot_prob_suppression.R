suppressMessages(library(tidyverse))
suppressMessages(library(patchwork))
suppressMessages(source('scripts/utils.R'))


pred_plot = function(dat, recs){
	p = ggplot() + 
		geom_rect(data=recs, aes(xmin=starts, xmax=ends, ymin=-2.5, ymax=102.5), fill='#eaeaea')
	for (s in c('susceptible', 'intermediate/high')){
		p = p + 
			geom_line(data=dat %>% filter(status == s),
				aes(x=round_mid_date, y=median*100, group=status, color=status),
				linewidth=2, color='white') + 
			geom_point(data=dat %>% filter(status == s),
				aes(x=round_mid_date, y=median*100, group=status, color=status),
				fill='white', color='white', shape=21, size=5) +
			geom_errorbar(data=dat %>% filter(status == s),
				aes(x=round_mid_date, ymin=`095lower`*100 -0.25, ymax=`095upper`*100 + 0.25),
				color='white', linewidth=2, width=0) +
			geom_errorbar(data=dat %>% filter(status == s),
				aes(x=round_mid_date, ymin=`050lower`*100 - 0.25, ymax=`050upper`*100 + 0.25),
				color='white', linewidth=3, width=0) +
			geom_errorbar(data=dat %>% filter(status == s),
				aes(x=round_mid_date, ymin=`095lower`*100, ymax=`095upper`*100, color=status),
				linewidth=1.5, width=0) +
			geom_errorbar(data=dat %>% filter(status == s),
				aes(x=round_mid_date, ymin=`050lower`*100, ymax=`050upper`*100, color=status),
				linewidth=2.5, width=0) +
		 	geom_line(data=dat %>% filter(status == s),
				aes(x=round_mid_date, y=median*100, group=status, color=status), 
		 		linewidth=1.5) + 
		 	geom_point(data=dat %>% filter(status == s),
				aes(x=round_mid_date, y=median*100, group=status, color=status),
		 		shape=21, fill='#eaeaea', size=3, stroke=1.5)
	}
	return(p + scale_y_continuous(limits=c(-2.5, 102.5), expand=c(0,0), name='posterior prob.\nof suppression (%)') + 
		 	xlab('date') +
		 	guides(color=guide_legend(position = "inside")) +
		 	gtheme +
		 	theme(legend.justification.inside = c(0, 1),
				legend.position.inside=c(0,1),
				panel.border = element_blank(),
				axis.line = element_line(colour = "#333333"),
				panel.grid.major=element_blank(),
				panel.grid.minor=element_blank(),
				legend.key = element_blank()))
}


rr_plot = function(rr, recs){
	return(
		ggplot(rr) +
			geom_hline(aes(yintercept=1), color='#848484', linetype='dashed') +
			#geom_rect(data=recs, aes(xmin=starts, xmax=ends, ymin=0, ymax=1), fill='#eaeaea') +
			geom_line(aes(x=round_mid_date, y=median),
				linewidth=2, color='white') + 
			geom_point(aes(x=round_mid_date, y=median),
				fill='white', color='white', shape=21, size=5) +
			geom_errorbar(
				aes(x=round_mid_date, ymin=`095lower`, ymax=`095upper`),
				color='white', linewidth=2, width=0) +
			geom_errorbar(
				aes(x=round_mid_date, ymin=`050lower`, ymax=`050upper`),
				color='white', linewidth=3, width=0) +
			geom_errorbar(
				aes(x=round_mid_date, ymin=`095lower`, ymax=`095upper`),
				linewidth=1.5, width=0, color='#333333') +
			geom_errorbar(
				aes(x=round_mid_date, ymin=`050lower`, ymax=`050upper`),
				linewidth=2.5, width=0, color='#333333') +
		 	geom_line(aes(x=round_mid_date, y=median), 
		 		linewidth=1.5, color='#333333') + 
		 	geom_point(aes(x=round_mid_date, y=median),
		 		shape=21, fill='#eaeaea', size=3, stroke=1.5, color='#333333') +
		 	ylim(0,2) +
		 	ylab('posterior prob. of\nsuppresison ratio') +
		 	xlab('date') +
		 	gtheme +
		 	theme(panel.border = element_blank(),
				axis.line = element_line(colour = "#333333"),
				panel.grid.major=element_blank(),
				panel.grid.minor=element_blank(),
				legend.key = element_blank()))
}


dodge = c('NA' = 0, 'susceptible' = -20, 'intermediate/high' = 20)


all_pred = read_tsv('models/overall_posterior_prob_suppression.tsv', show_col_types=FALSE) %>%
		filter(grepl('round_mu_pred\\[', name) & status == "susceptible") %>%
		mutate(round_mid_date = as.Date(round_mid_date))

nnrti = read_tsv('models/nnrti_posterior_prob_suppression.tsv', show_col_types=FALSE) %>%
	filter(grepl('round_mu_pred\\[', name) | grepl('within_round_rr\\[', name))  %>%
		mutate(round_mid_date = if_else(
			is.na(status), 
				as.Date(round_mid_date),
				as.Date(round_mid_date)+ dodge[as.character(status)]))

nrti = read_tsv('models/nrti_posterior_prob_suppression.tsv', show_col_types=FALSE) %>%
	filter(grepl('round_mu_pred\\[', name) | grepl('within_round_rr\\[', name))  %>%
		mutate(round_mid_date = if_else(
			is.na(status), 
				as.Date(round_mid_date),
				as.Date(round_mid_date)+ dodge[as.character(status)]))


d = read_tsv('data/rakai_drug_resistance_categorized_R15_R20.tsv', show_col_types=FALSE)  %>%
	select(round, round_min_date, round_mid_date, round_max_date) %>%
	unique()
recs = as_tibble(d %>% slice(seq(nrow(d),1,-2)) %>%
	select(round_min_date, round_max_date) %>%
	rename(c('starts' = 'round_min_date',
		'ends' = 'round_max_date')) %>%
	mutate(starts = as.Date(starts), ends=as.Date(ends)))

pA = pred_plot(nnrti %>% filter(grepl('round_mu_pred\\[', name)), recs) +
	scale_color_manual(
		values=c('susceptible'='#333333', 'intermediate/high'=colors[['nnrti']]),
		labels=c('susceptible'='susceptible', 'intermediate/high' = "nnrti"),
		name=NULL) 

pB = rr_plot(nnrti %>% filter(grepl('within_round_rr\\[', name)), recs) + 
	ggtitle('nnrti')

pC = pred_plot(nrti %>% filter(grepl('round_mu_pred\\[', name)), recs) +
	scale_color_manual(
		values=c('susceptible'='#333333', 'intermediate/high'=colors[['nrti']]),
		labels=c('susceptible'='susceptible', 'intermediate/high' = "nrti"),
		name=NULL)

pD = rr_plot(nrti %>% filter(grepl('within_round_rr\\[', name)), recs) + 
	ggtitle('nrti')
 		
p = (pA + pB) / (pC + pD) + 
	plot_annotation(
		tag_levels = list(c('A', 'B', 'C', 'D'))) &
	theme(
		plot.tag = element_text(color='#333333', size=16, face='bold'))

ggsave('figures/pdf/p_suppression_among_viremic.pdf', p, width=6.4*1.75, height=4.2*1.75)
ggsave('figures/png/p_suppression_among_viremic.png', p, width=6.4*1.75, height=4.2*1.75)
ggsave('figures/poster_p_suppression_among_viremic.pdf', pC +xlab(NULL), width=4.2, height=3.6)



pA = pred_plot(all_pred %>% mutate(status = "susceptible"), recs) +
	scale_color_manual(
		values=c('susceptible'='#333333')) + 
	guides(color="none")
ggsave('figures/pdf/p_suppression_among_viremic_all.pdf', pA, width=6.4, height=4.2)
ggsave('figures/png/p_suppression_among_viremic_all.png', pA, width=6.4, height=4.2)

