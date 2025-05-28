suppressMessages(library(tidyverse))
suppressMessages(library(tidytree))
suppressMessages(library(treeio))
suppressMessages(require(ggtree))
suppressMessages(require(ape))
suppressMessages(library(shadowtext))
suppressMessages(library(patchwork))
suppressMessages(source('scripts/utils.R'))

get_lower_tri = function(x){
	x[upper.tri(x)] = NA
	return(x)}


muts = bind_rows(
	read_tsv('data/rakai_drug_resistance_mut_R15_R20.tsv', show_col_types=FALSE),
	read_tsv('data/other_rakai_drug_resistance_mut.tsv', show_col_types=FALSE))

ins153y_freq = muts %>% 
	filter(mut == 'inS153Y') %>% 
	select(study_id, int_date, mut, freq) %>% 
	unique()

# manually bin it
p = ggplot(ins153y_freq, aes(x=freq*100)) +
	geom_vline(aes(xintercept=5), linetype='dashed', color='#848484') +
	geom_histogram(color='#333333', fill=colors['insti'], binwidth = 1, center=0.5) +
	xlab('inS153Y frequency') + 
	ylab('participants') +
	xlim(0,100) +
	scale_y_continuous(expand=c(0,0), limits=c(0,15)) +
	gtheme +
	theme(panel.border = element_blank(),
		axis.line = element_line(colour = "#333333"))

ggsave('figures/pdf/ins153y_freq.pdf', p, width=6.4, height=4.8)
ggsave('figures/png/ins153y_freq.png', p, width=6.4, height=4.8)

