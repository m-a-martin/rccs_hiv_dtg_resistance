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

dat = bind_rows(
	read_tsv('data/rakai_drug_resistance_categorized_R15_R20.tsv', show_col_types=FALSE),
	read_tsv('data/other_rakai_drug_resistance_categorized.tsv', show_col_types=FALSE)) %>%
	filter(viremic & !is.na(insti) & !is.na(nnrti) & !is.na(nrti) & !is.na(pi)) %>%
	select(study_id, round, sex, int_date, Mut, numeric_copies, viremic, subtype_bestref)

# get just viremic dat after start of R20
r20_dat = dat %>% filter(int_date >= "2021-02-08") %>%
	mutate(
		inS153Y = grepl('inS153Y', Mut),
		seq=paste(study_id, str_split(as.character(int_date), " ", simplify=TRUE)[,1], sep='_'))

# read in dists
dists = read_csv('data/2024-10-02_pangea2_aln_pol_mask_dist.csv', show_col_types=FALSE) %>%
	filter(seq1 != seq2)
dists = dists %>%
	inner_join(r20_dat %>% select(seq, inS153Y), by=c('seq1'='seq')) %>%
	inner_join(r20_dat %>% select(seq, inS153Y), by=c('seq2'='seq'))

# 99th percentile of non-inS153Y dists
thresh = (dists %>% filter(inS153Y.x == FALSE & inS153Y.y == FALSE) %>% summarise(x = quantile(d, 0.01)))$x
print(thresh)

dists %>% filter(inS153Y.x == TRUE & inS153Y.y == TRUE) %>% arrange(d)

# plot inS153Y dists
p1 = ggplot(dists %>% filter(inS153Y.x == TRUE & inS153Y.y == TRUE), aes(x=d)) +
	geom_vline(aes(xintercept=thresh), linetype='dashed', color='#848484') +
	geom_histogram(color='#333333', fill=colors['insti'], binwidth = 0.001) +
	xlab('genetic distance (substitutions / site)') + 
	ylab('pairwise comparisons') +
	scale_y_continuous(expand=c(0,0)) +
	gtheme +
	theme(panel.border = element_blank(),
		axis.line = element_line(colour = "#333333"))

# plot A1 tree
a1 = read.newick('data/2024-10-02_pangea2_aln_pol_mask_A1_inS153Y.fasta.treefile')
a1 = drop.tip(a1, "K03455.1")
#a1 = as.treedata(a1)
g = '#333333'
ymax1 = 0.0905
p2 = ggplot() +
		geom_tile(aes(
			x=seq(ymax1/2,ymax1,0.0005), 
			y=rep(37.5,length(seq(ymax1/2,ymax1,0.0005))),
			alpha=seq(0,1,1/(length(seq(ymax1/2,ymax1,0.0005))-1))),
			fill=colors['insti'],
			height=5) +
		geom_tree(data=a1, aes(x=x,y=y), color='white', linewidth=2) +
		geom_tree(data=a1, aes(x=x,y=y), color=g, linewidth=1) +
		geom_tippoint(data=a1,aes(x=x,y=if_else(grepl('True', label), y, NA)),
			color="white", fill="white",  shape=21, stroke=1, size=4) +
		geom_tippoint(data=a1, aes(x=x,y=if_else(grepl('True', label), y, NA)),
			fill=colors['insti'], color=g, shape=21, stroke=1, size=3) +
		geom_shadowtext(
	    	data = tibble(x=c(0.0125), y=c(nrow(as_tibble(a1) %>% filter(!is.na(label)))*0.9),
    			label=c("inS153Y")),
	   		aes(x, y,  label = label),
	   		color=colors['insti'],
	   		size = 6,
	    	bg.color = g, bg.r=0.1) +
		ggtitle('subtype A1') +
		scale_alpha_continuous(range=c(0,0.75), guide="none") +
		scale_y_continuous(breaks=NULL, name=NULL) +
		xlab('substitutions/site') +
		gtheme +
		theme(
			plot.title=element_text(size=20),
			panel.grid.major.y = element_blank(),
			panel.grid.minor.y = element_blank(),
			axis.line.x = element_line(color='#333333'),
			panel.border = element_blank(),
			plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 

# d tree
d = read.newick('data/2024-10-02_pangea2_aln_pol_mask_D_inS153Y.fasta.treefile')
d = drop.tip(d, "K03455.1")
#a1 = as.treedata(a1)
g = '#333333'
ymax2 = 0.085
p3 = ggplot() +
		geom_tile(aes(
			x=seq(ymax2/2,ymax2,0.0005), 
			y=rep(39,length(seq(ymax2/2,ymax2,0.0005))),
			alpha=seq(0,1,1/(length(seq(ymax2/2,ymax2,0.0005))-1))), 
			height=4, fill=colors['insti']) +
		geom_tree(data=d, aes(x=x,y=y), color='white', linewidth=2) +
		geom_tree(data=d, aes(x=x,y=y), color=g, linewidth=1) +
		geom_tippoint(data=d,aes(x=x,y=if_else(grepl('True', label), y, NA)),
			color="white", fill="white", size=4, shape=21, stroke=1) +
		geom_tippoint(data=d, aes(x=x,y=if_else(grepl('True', label), y, NA)),
			fill=colors['insti'], color=g, size=3, shape=21, stroke=1) +
		ggtitle('subtype D') +
		scale_alpha_continuous(range=c(0,0.75), guide="none") +
		scale_y_continuous(breaks=NULL, name=NULL) +
		xlab('substitutions/site') +
		gtheme +
		theme(
			plot.title=element_text(size=20),
			panel.grid.major.y = element_blank(),
			panel.grid.minor.y = element_blank(),
			axis.line.x = element_line(color='#333333'),
			panel.border = element_blank(),
			plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 


p = (p1 + p2 + p3) + 
	plot_annotation(
		tag_levels = list(c('A', 'B', 'C'))) & 
	theme(
		plot.tag = element_text(color='#333333', size=16, face='bold'))

ggsave('figures/pdf/ins153y_pairwise_dists.pdf', p, width=14, height=4.8)
ggsave('figures/png/ins153y_pairwise_dists.png', p, width=14, height=4.8)
ggsave('figures/poster_ins153y_pairwise_dists.pdf', p1, width=4.2, height=3.6)


