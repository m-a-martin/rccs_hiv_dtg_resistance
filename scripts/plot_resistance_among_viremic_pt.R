suppressMessages(library(tidyverse))
suppressMessages(library(geepack))
suppressMessages(library(emmeans))
suppressMessages(library(patchwork))
suppressMessages(library(cowplot))
suppressMessages(source('scripts/utils.R'))
 
#dodge = c('any'=-40, 'insti' = -20, 'nnrti'=0, "nrti" = 20, "pi" = 40)
dodge = c('insti' = -30, 'nnrti'=-10, "nrti" = 10, "pi" = 30)

## PANEL A ##
pred = bind_rows(
	#read_tsv('models/any_among_viremic_pt_pred.tsv', show_col_types=FALSE),
	read_tsv('models/insti_among_viremic_pt_pred.tsv', show_col_types=FALSE),
	read_tsv('models/nnrti_among_viremic_pt_pred.tsv', show_col_types=FALSE),
	read_tsv('models/nrti_among_viremic_pt_pred.tsv', show_col_types=FALSE),
	read_tsv('models/pi_among_viremic_pt_pred.tsv', show_col_types=FALSE)) %>%
	mutate(label = str_split(label, '_', simplify=TRUE)[,1]) %>%
	mutate(label = ordered(label, levels=unique(label)),
		round_mid_date = as.Date(round_mid_date)+ dodge[label]) %>%
	mutate(lwr = if_else(obs == 0, NA, lwr),
		upr = if_else(obs == 0, NA, upr))


#	group_by(round) %>%
#	summarise(
#		round_min_date = as.Date(quantile(int_date, 0.01, na.rm=TRUE)),
#		round_mid_date = as.Date(median(int_date, na.rm=TRUE)),
#		round_max_date = as.Date(quantile(int_date, 0.99, na.rm=TRUE)),
#		.groups='drop') %>%
#	merge(pred %>% select(round) %>% unique(),
#		how='inner', on='round')

d = read_tsv('data/rakai_drug_resistance_categorized_R15_R20.tsv', show_col_types=FALSE)  %>%
	select(round, round_min_date, round_mid_date, round_max_date) %>%
	unique()
recs = as_tibble(d %>% slice(seq(nrow(d),1,-2)) %>%
	select(round_min_date, round_max_date) %>%
	rename(c('starts' = 'round_min_date',
		'ends' = 'round_max_date')) %>%
	mutate(starts = as.Date(starts), ends=as.Date(ends)))

ymax=5*ceiling(max(pred$upr, na.rm=TRUE)*100/5)
pA = ggplot(pred) + 
	geom_rect(data=recs, aes(xmin=starts, xmax=ends, ymin=-ymax*0.03, ymax=ymax), fill='#eaeaea') +
	geom_line(
		aes(x=round_mid_date, y=fit*100, group=label),
		linewidth=2, color='white') + 
	geom_point(
		aes(x=round_mid_date, y=fit*100, group=label, shape=label),
		fill='white', color='white', size=5) +
	geom_errorbar(
		aes(x=round_mid_date, ymin=lwr*100, ymax=upr*100, group=label),
		color='white', linewidth=2, width=0) +
	geom_errorbar(
		aes(x=round_mid_date, ymin=lwr*100, ymax=upr*100, group=label, color=label),
		width=0, linewidth=1.5) +
	geom_line(
		aes(x=round_mid_date, y=fit*100, group=label, color=label),
		linewidth=1.5) + 
	geom_point(
		aes(x=round_mid_date, y=fit*100, group=label, color=label, shape=label),
		fill='#eaeaea', size=3, stroke=1.5) +
	scale_y_continuous(expand=c(0,0), limits=c(-ymax*0.03, ymax)) +
	ylab('prevalence (%)') +
	xlab('date') + 
	scale_shape_manual(values=shapes, name=NULL) +
	scale_color_manual(values=colors, name=NULL) +
	guides(
		shape=guide_legend(position="inside"),
		color=guide_legend(position="inside")) +
	gtheme +
	theme(
		legend.justification.inside = c(0, 1),
		legend.position.inside=c(0,1),
		panel.border = element_blank(),
		axis.line = element_line(colour = "#333333"),
		panel.grid.major=element_blank(),
		panel.grid.minor=element_blank(),
		legend.key = element_blank()) 


### PANEL B ###
files = Sys.glob(file.path('models/mutations', "*-pt-*pred.tsv"))
mut_pred = do.call('rbind',lapply(files,read_tsv)) %>% 
	filter(round == max(round)) %>%
	mutate(mut = str_split(label, "-", simplify=TRUE)[,3]) %>% 
	arrange(-fit) %>%
	left_join(
		read_tsv('config/mut_class.tsv', show_col_types=FALSE),
		by='mut') %>%
	mutate(mut = ordered(mut, levels=mut))

pB = ggplot(mut_pred) + 
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


### PANEL C #####
# todo fix correlatino structure
dodge = c('rtK103N' = -10, 'rtE138A' = 10)
files = c('models/mutations/viremic-pt-rtK103N_pred.tsv',
	'models/mutations/viremic-pt-rtE138A_pred.tsv')
long_mut_pred = do.call('rbind',lapply(files,function(x){read_tsv(x, show_col_types=FALSE)})) %>% 
	mutate(mut = str_split(label, "-", simplify=TRUE)[,3]) %>%
	arrange(-fit) %>%
	left_join(
		read_tsv('config/mut_class.tsv', show_col_types=FALSE),
		by='mut') %>%
	mutate(
		label = ordered(mut, levels=unique(mut)),
		round_mid_date = as.Date(round_mid_date)+ dodge[label])


ymax=5*ceiling(max(pred$upr, na.rm=TRUE)*100/5)
pC = ggplot(long_mut_pred) + 
	geom_rect(data=recs, aes(xmin=starts, xmax=ends, ymin=-ymax*0.03, ymax=ymax), fill='#eaeaea') +
	geom_line(
		aes(x=round_mid_date, y=fit*100, group=label),
		linewidth=2, color='white') + 
	geom_point(
		aes(x=round_mid_date, y=fit*100, group=label, shape=label),
		fill='white', color='white', size=5) +
	geom_errorbar(
		aes(x=round_mid_date, ymin=lwr*100, ymax=upr*100, group=label),
		color='white', linewidth=2, width=0) +
	geom_errorbar(
		aes(x=round_mid_date, ymin=lwr*100, ymax=upr*100, group=label, color=label),
		width=0, linewidth=1.5) +
	geom_line(
		aes(x=round_mid_date, y=fit*100, group=label, color=label),
		linewidth=1.5) + 
	geom_point(
		aes(x=round_mid_date, y=fit*100, group=label, color=label, shape=label),
		fill='#eaeaea', size=3, stroke=1.5) +
	scale_y_continuous(expand=c(0,0), limits=c(-ymax*0.03, ymax)) +
	ylab('prevalence (%)') +
	xlab('date') + 
	scale_shape_manual(values=shapes, name=NULL) +
	scale_color_manual(values=colors, name=NULL) +
	guides(
		shape=guide_legend(position="inside"),
		color=guide_legend(position="inside")) +
	gtheme +
	theme(
		legend.justification.inside = c(0, 1),
		legend.position.inside=c(0,1),
		panel.border = element_blank(),
		axis.line = element_line(colour = "#333333"),
		panel.grid.major=element_blank(),
		panel.grid.minor=element_blank(),
		legend.key = element_blank()) 


p1 = (pA + 	theme(plot.tag = element_text(color='#333333', size=16, face='bold'))) / (pC + 
	theme(plot.tag = element_text(color='#333333', size=16, face='bold')))
p1 = p1 + 
	plot_annotation(tag_levels = list(c('A', 'C'))) 

p2 = (plot_spacer() / (pB +
	theme(plot.tag = element_text(color='#333333', size=16, face='bold'))) / plot_spacer()) + plot_layout(heights = c(0.5, 1, 0.5)) + 
	plot_annotation(tag_levels = list(c('B'))) 

title <- ggdraw() + draw_label("pre-treatment viremic PLHIV", size=22, color='#333333')
p = plot_grid(p1, p2, rel_widths=c(1.5,1))
p = plot_grid(title, p, ncol=1, rel_heights=c(0.05, 1))

ggsave('figures/pdf/resistance_among_viremic_pt.pdf', width=6.4*1.9, height=4.8*1.9)
ggsave('figures/png/resistance_among_viremic_pt.png', width=6.4*1.9, height=4.8*1.9)


p = (pA + xlab(NULL) + theme(legend.position.inside=c(0,1), legend.justification.inside=c(0,1))) + 
		plot_spacer()  + pB + plot_layout(widths = c(1.5, 0.15, 1))+ 
	plot_annotation(
		title = 'pre-treatment PLHIV')&
	theme(
		plot.tag = element_text(color='#333333', size=16, face='bold'),
		plot.title = element_text(color='#333333', size=22, hjust = 0.5))
ggsave('figures/poster_resistance_among_viremic_pt.pdf', p, width=10.8, height=4.2)

